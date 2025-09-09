import re

rule fetch_ecoligenes:
  output:
    source = "%s/custom/ecoligenes.fasta" %database_path
  conda:
    "../envs/fetch.yaml"
  log:
    stdout = "Logs/Databases/setup_ecoligenes_ecoligenes.log"
  message:
    "[fetch_ecoligenes]: Downloading custom database ecoligenes"
  shell:
    """
    mkdir -p $(dirname {output.source})
    cmd="curl https://raw.githubusercontent.com/ssi-dk/ecoli_fbi/refs/heads/main/db/ecoligenes/ecoligenes.fsa -o {output.source}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1

    """


rule samtools_faidx_directory:
    input:
        fasta_dir = lambda wc: "%s/%s" %(database_path, db_flat_to_path(wc.db_name_flat))
    output:
        flag = "%s/{db_name_flat}.fa_idx.done" %database_path
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[samtools_faidx_directory]: Indexing all FASTA in {input.fasta_dir}"
    shell:
        """
        mkdir -p {input.fasta_dir}
    
        for fasta in $(find {input.fasta_dir} -maxdepth 1 -type f  -name "*.f*a*"); do
            prefix="${{fasta%.*}}"
            echo "Indexing $fasta -> $prefix"
            cmd="samtools faidx $fasta"
            echo "Executing command:\n$cmd\n" >> {output.flag} 
            eval $cmd >> {output.flag} 2>&1
        done

        date -I >> {output.flag}
        """


rule setup_Ecoli_alignment_db:
    conda:
        config["analysis_settings"]["Escherichia_coli_db"]["yaml"]
    params:
        db_prefix = 'ecoligenes'
    output:
        database = directory("%s/%s" % (database_path, config["analysis_settings"]["Escherichia_coli_db"]["database"]))
    log:
        stdout = 'Logs/Databases/setup_Ecoli_alignment_db.log'
    message:
        "[setup_Ecoli_alignment_db]: Setting up Ecoli alignment database"
    shell:
        """
        mkdir -p {output.database}
        cmd="curl https://raw.githubusercontent.com/ssi-dk/ecoli_fbi/refs/heads/main/db/ecoligenes/ecoligenes.fsa -o {output.database}/{params.db_prefix}.fsa"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        date -I > {output.database}/creation.date
        """


rule setup_CdiffToxin:
    conda:
        config["analysis_settings"]["Clostridioides_difficile_db"]["yaml"]
    output:
        database = directory("%s/%s/Toxin" % (database_path, config["analysis_settings"]["Clostridioides_difficile_db"]["database"]))
    params:
        accession_loci = "AM180355.1:tcdA,tcdB,tcdC AF271719.1:cdtA,cdtB",
        db_toxin = "Cdiff_Toxin"
    log:
        stdout = 'Logs/Databases/setup_CdiffToxin.log'
    message:
        "[setup_CdiffToxin]: Setting up C. difficile toxin database"
    shell:
        """
        set -euo pipefail
        mkdir -p {output.database}
        > {output.database}/{params.db_toxin}.txt
        > {output.database}/{params.db_toxin}.bed6
        > {output.database}/{params.db_toxin}.fasta

        for item in {params.accession_loci}; do
            acc=$(echo $item | cut -d':' -f1)
            loci=$(echo $item | cut -d':' -f2 | tr ',' ' ')

            cmd="python workflow/scripts/genbank_fetcher.py \
                -a $acc --locus $loci \
                -o {output.database}/{params.db_toxin}.txt \
                --bed {output.database}/{params.db_toxin}.bed6 \
                --fasta {output.database}/{params.db_toxin}.fasta \
                --merge 500 --append"

            echo -e \"\nExecuting command:\n$cmd\n\" >> {log.stdout}
            eval $cmd >> {log.stdout} 2>&1
        done
        date -I > {output.database}/creation.date
        """

######### DATABASE REPEAT PATTERNS FOR CDIFF

species_key = config["analysis_settings"]["Clostridioides_difficile_db"]["species_config"]
repeat_list = species_configs[species_key]["analyses_to_run"]["Repeat_identifier"]["repeats"].split()
combo_list = species_configs[species_key]["analyses_to_run"]["Repeat_identifier"]["combos"].split()
TR_repeat_sequences_str = " ".join(repeat_list)
TR_repeat_types_str = " ".join(repeat_list + combo_list)

trst_output_repeat_fa = expand(
    "%s/%s/TRST/{repeat}_repeat_sequences.fa" % (database_path, config["analysis_settings"]["Clostridioides_difficile_db"]["database"]),
    repeat=repeat_list
)

rule setup_CdiffTRST:
    conda:
        config["analysis_settings"]["Clostridioides_difficile_db"]["yaml"]
    output:
        database = directory("%s/%s/TRST" % (database_path, config["analysis_settings"]["Clostridioides_difficile_db"]["database"])),
        repeat_fa = trst_output_repeat_fa
    params:
        TR_repeat_sequences = TR_repeat_sequences_str,
        TR_repeat_types = TR_repeat_types_str,
    log:
        stdout = 'Logs/Databases/setup_CdiffTRST.log'
    message:
        "[setup_CdiffTRST]: Downloading TRST repeat sequences and types"
    shell:
        """
        set -euo pipefail
        mkdir -p {output.database}

        # TR-type repeat sequence downloads
        for TR in {params.TR_repeat_sequences}; do
            cmd="curl -fSL https://raw.githubusercontent.com/ssi-dk/cdiff_fbi/refs/heads/raah_dev/db/TRST/${{TR}}_repeat_sequences.fa \
                -o {output.database}/${{TR}}_repeat_sequences.fa"
            echo "Executing command:\n$cmd\n" >> {log.stdout}
            eval $cmd >> {log.stdout} 2>&1
        done

        # TR-type repeat types downloads
        for TR in {params.TR_repeat_types}; do
            cmd="curl -fSL https://raw.githubusercontent.com/ssi-dk/cdiff_fbi/refs/heads/raah_dev/db/TRST/${{TR}}_repeat_types.txt \
                -o {output.database}/${{TR}}_repeat_types.txt"
            echo "Executing command:\n$cmd\n" >> {log.stdout}
            eval $cmd >> {log.stdout} 2>&1
        done
        """



rule setup_all_databases:
    input:
        rules.setup_Ecoli_alignment_db.output.database,
        rules.setup_CdiffToxin.output.database,
        rules.setup_PlasmidFinder.output.database,
        rules.setup_ResFinder.output.database,
        rules.setup_PointFinder.output.database,
        rules.setup_DisinFinder.output.database,
        rules.setup_VirulenceFinder.output.database,
        rules.setup_AMRFinder.output.database,
        rules.setup_CdiffTRST.output.database
