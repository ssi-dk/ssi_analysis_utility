#!/usr/bin/env python3
"""
spaTyper v2.0

Reimplementation of the spaTyper tool from CGE: https://bitbucket.org/genomicepidemiology/spatyper
to predicts the S.Aureus spa type from genome sequences.

Contributors:
- Simone Scrima
"""

# =========================
# Imports
# =========================
import os
import glob
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import pandas as pd
from Bio import SeqIO
import argparse
import sys
from git import Repo


# ============================
# Logging class and functions
# ============================


class Logging_ArgParser(argparse.ArgumentParser):
    def error(self, message):
        logging.error(message)
        self.print_usage(sys.stderr)
        sys.exit(2)


def setup_logging(logfile):
    logging.basicConfig(
        filename=logfile,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    

# =========================
# Utility functions
# =========================

def batch_iterator(iterator, 
                   batch_size):
    """
    Yield batches of records from an iterator.

    Parameters
    ----------
    iterator : iterator
        Iterator over records (e.g. SeqIO.parse).
    batch_size : int
        Number of records per batch.

    Yields
    ------
    list
        List of records.
    """
    batch = []
    for record in iterator:
        batch.append(record)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def clean_tmp():
    """Remove temporary FASTA and BLAST output files."""
    for tmp in glob.glob("tmp_spaseq_*.fasta"):
        os.remove(tmp)
    for tmp in glob.glob("blastn_out_tmp_spaseq*.tab"):
        os.remove(tmp)
    os.remove("combined.tab")


# =========================
# spaTyper class
# =========================

class SpaTyper:
    """
    spaTyper pipeline class.
    """

    def __init__(self, 
                 assembly, 
                 spatype_db, 
                 blast_db, 
                 ncore):
        
        self.assembly = assembly
        self.spatype_db = spatype_db
        self.blast_db = blast_db
        self.ncore = ncore
    
        self.combined_output = "combined.tab"
        self.saco_dict = SeqIO.index(self.assembly, "fasta")
        

    # ---------------------
    # FASTA handling
    # ---------------------

    def split_fasta(self):
        """Split spa FASTA file into chunks for parallel BLAST."""
        
        self.spa_sequences = os.path.join(self.spatype_db, "spa_sequences.fna")
        record_dict = SeqIO.index(self.spa_sequences, "fasta")
        batch_size = max(1, int(len(record_dict) / self.ncore))

        record_iter = SeqIO.parse(self.spa_sequences, "fasta")

        for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
            out = f"tmp_spaseq_{i + 1}.fasta"
            with open(out, "w") as handle:
                SeqIO.write(batch, handle, "fasta")



    # =========================
    # Databases validations
    # =========================

    def validate_blast_db(self):
        """
        Validate the BLAST database directory.
        If no valid nucleotide BLAST DB is found, create one from the assembly.
        """

        # Minimum required files for a nucleotide BLAST DB
        REQUIRED_BLAST_EXT = {".ndb",
                              ".nin",
                              ".not",
                              ".ntf",
                              ".nhr",
                              ".njs",
                              ".nsq",
                              ".nto"}

        # Directory exists
        os.makedirs(self.blast_db, exist_ok=True)

        files = os.listdir(self.blast_db)

        # Extensions found
        prefix_map = {}

        for name in files:
            base, ext = os.path.splitext(name)
            if ext in REQUIRED_BLAST_EXT:
                prefix_map.setdefault(base, set()).add(ext)

        # Check for any valid DB using sets
        for prefix, exts in prefix_map.items():
            if REQUIRED_BLAST_EXT.issubset(exts):
                detected_prefix = os.path.join(self.blast_db, prefix)
                logging.info("Detected existing BLAST DB: %s", detected_prefix)
                self.blast_db = detected_prefix
                return self.blast_db

        # No valid DB found -> create one
        logging.warning(
            "No valid BLAST nucleotide database found in %s. Creating from assembly.",
            self.blast_db)

        out_prefix = os.path.join(self.blast_db,
                                  os.path.basename(os.path.splitext(self.assembly)[0]))

        self._run_makeblastdb(self.assembly, out_prefix)

        # Re-check if all is good
        missing = []
        for ext in REQUIRED_BLAST_EXT:
            if not os.path.isfile(f"{out_prefix}{ext}"):
                missing.append(f"{out_prefix}{ext}")

        if missing:
            logging.error(
                "BLAST DB creation failed. Missing files:\n  " +
                "\n  ".join(missing)
            )
            raise RuntimeError("makeblastdb failed")

        logging.info("BLAST DB successfully created: %s", out_prefix)
        self.blast_db = out_prefix
        return self.blast_db


    def validate_spatype_db(self):
        """
        Validate the spatype database.
        spatype_db is the value passed to -d /--spatype_db
        """

        REQUIRED_SPATYPE_FILES = ("spa_types.txt",
                                  "spa_sequences.fna")
        missing_spatype = []
        
        for file in REQUIRED_SPATYPE_FILES:
            file = os.path.join(self.spatype_db, file)
            if not os.path.isfile(f"{file}"):
                missing_spatype.append(f"{file}")
        
        
        # If missing report it and download it in the current spatype_db directory path
        if missing_spatype:
            logging.error(
                "Invalid spatype database.\n"
                f"Missing files:\n  " + "\n  ".join(missing_spatype)) 
            Repo.clone_from("https://git@bitbucket.org/genomicepidemiology/spatyper_db.git", 
                            self.spatype_db)
            for file in REQUIRED_SPATYPE_FILES:
                if not os.path.isfile(f"{self.spatype_db}{file}"):
                    raise RuntimeError(
                        f"Spatyper clone failed: {self.spatype_db}{file} not downloaded")
            logging.info("Spatyper database successfully downloaded: %s", self.spatype_db)        
            return self.spatype_db
        
        return self.spatype_db 


    # ---------------------
    # BLAST
    # ---------------------

    def _run_makeblastdb(self,
                         fasta,
                         db_path):
        """Run makeblast_db  command on a fasta file"""

        cmd = [
            "makeblastdb",
            "-in", fasta,
            "-dbtype", "nucl",
            "-out", db_path
        ]    

        subprocess.run(cmd, check=True)

    def _run_blast(self, 
                   query_file):
        """Run BLASTn on a single FASTA chunk."""
        
        base = os.path.splitext(os.path.basename(query_file))[0]
        output = f"blastn_out_{base}.tab"

        cmd = [
            "blastn",
            "-query", query_file,
            "-db", self.blast_db,
            "-out", output,
            "-outfmt", "6",
            "-dust", "no",
            "-evalue", "0.0001",
            "-num_alignments", "10000"
        ]

        subprocess.run(cmd, check=True)
        return output

    def run_blast_parallel(self):
        """Run BLASTn in parallel and merge outputs."""
        input_files = glob.glob("tmp_spaseq_*.fasta")
        outputs = []

        with ThreadPoolExecutor(max_workers=self.ncore) as executor:
            
            futures = []
            for f in input_files:
                futures.append(executor.submit(self._run_blast, f))

            for future in as_completed(futures):
                outputs.append(future.result())
                logging.info(f"Finished: {future.result()}")

        with open(self.combined_output, "w") as out:
            for f in sorted(outputs):
                with open(f) as fh:
                    out.write(fh.read())
            logging.info(f"Combined BLAST output written to: {self.combined_output}")

        
        

    # ---------------------
    # Parsing & filtering
    # ---------------------

    def load_and_filter_hits(self):
        """
        Load BLAST output and apply quality filters.

        Returns
        -------
        pandas.DataFrame
        """
        columns = [
            "qseqid", 
            "sseqid", 
            "pident", 
            "length", 
            "mismatch", 
            "gapopen",
            "qstart", 
            "qend", 
            "sstart", 
            "send", 
            "evalue", 
            "bitscore"
        ]

        df = pd.read_csv(self.combined_output, 
                         sep="\t", 
                         names=columns)

        return df[
            (df["pident"] >= 100.0) &
            (df["length"] >= 20) &
            ((df["qstart"] == 1) | (df["qend"] == 1))
        ]

    # ---------------------
    # spa repeat logic
    # ---------------------

    def _load_repeats(self):
        """Load spa-type repeat database."""
        d_repeats = {}
        with open(os.path.join(self.spatype_db, "spa_types.txt")) as f:
            for line in f:
                type, repeat = line.strip().split(",")
                d_repeats[type] = repeat
        return d_repeats

    def match_spa_ends(self, 
                       filtered_df):
        """Validate spa ends using conserved flanking motifs."""
        
        results = []
        for _, row in filtered_df.iterrows():

            spatype = row["qseqid"].removeprefix("spatype_")
            contig = row["sseqid"]

            if contig not in self.saco_dict:
                continue

            qpos1, qpos2 = row["sstart"], row["send"]
            contig_seq = self.saco_dict[contig]

            
            if qpos1 < qpos2:
                orien = "plus" 
                pos1 = qpos1 - 10
                pos2 = qpos2 + 29
                end5 = contig_seq[pos2 - 10:pos2]
                end3 = contig_seq[pos1:qpos1 - 1]
            else:
                orien = "minus"
                pos1 = qpos2 - 30
                pos2 = qpos1 + 9
                end5 = contig_seq[pos1:pos1 + 10].reverse_complement()
                end3 = contig_seq[pos1:pos2].reverse_complement()

            if (
                re.search(r"TA[CT]ATGTCGT", end5.format("fasta")) and
                re.search(r"[GA]CA[AC]CAAAA", end3.format("fasta"))
            ):
                results.append([spatype,
                                self._load_repeats()[spatype],
                                contig,
                                pos1,
                                pos2,
                                orien])
        
        self.results_df = pd.DataFrame(results,
                                       columns=["spaType",
                                                "Repeats",
                                                "Contig",
                                                "Start",
                                                "End",
                                                "Orientation"])
        return self.results_df
    



# =========================
# Argument parser
# =========================

def parse_args():

    """Parse command-line arguments."""
    parser = Logging_ArgParser(
        description="spaTyper v2.0 — spa typing using BLAST and motif validation"
    )

    parser.add_argument(
        "-a", "--assembly",
        required=True,
        help="Assembly FASTA file"
    )

    parser.add_argument(
        "-d", "--spatype_db",
        required=True,
        help="Path to spaTyper database"
    )

    parser.add_argument(
        "-b", "--blast_db",
        default="seq_db",
        help="BLAST nucleotide database name"
    )

    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of parallel BLAST threads (default: 8)"
    )

    parser.add_argument(
        "-l", "--log",
        type=str,
        default="spaTyper.log",
        help="Output name"
    )

    parser.add_argument(
        "-o", "--out",
        type=str,
        default="spaTyper_results.tsv",
        help="Output name"
    )

    args = parser.parse_args()

    if not os.path.isfile(args.assembly):
        parser.error(f"Assembly file not found: {args.assembly}")

    logging.info("Command-line arguments: %s", vars(args))
    return args


def main():
    
    args = parse_args()
    setup_logging(args.log)

    spatyper = SpaTyper(assembly=args.assembly,
                        spatype_db=args.spatype_db,
                        blast_db=args.blast_db,
                        ncore=args.threads)



    spatyper.validate_blast_db()
    spatyper.validate_spatype_db()
    spatyper.split_fasta()
    spatyper.run_blast_parallel()
    
    filtered_hits = spatyper.load_and_filter_hits()
    results = spatyper.match_spa_ends(filtered_hits)
    results.to_csv(args.out,
                   sep="\t")
    clean_tmp()


if __name__ == "__main__":
    main()
