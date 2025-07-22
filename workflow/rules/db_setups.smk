import re

##### INDEXING #####
def db_flat_to_path(name):
    # Replace the first underscore after _db with a slash
    return re.sub(r'(_db)_(?=[^_]+$)', r'\1/', name)
    
rule kma_index_directory:
    output:
        flag = "Logs/Databases/{db_name_flat}.kma_index.done"
    input:
        fasta_dir = lambda wc: f"resources/{db_flat_to_path(wc.db_name_flat)}"
    conda:
        config["analysis_settings"]["KMA"]["yaml"]
    message:
        "[kma_index_directory]: Indexing all FASTA in {input.fasta_dir}"
    shell:
        r"""
        mkdir -p {input.fasta_dir}

        for fasta in $(find {input.fasta_dir} -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fsa" -o -name "*.fasta" \)); do
            prefix="${{fasta%.*}}"
            echo "Indexing $fasta -> $prefix"
            kma index -i "$fasta" -o "$prefix"
        done

        touch {output.flag}
        """

rule samtools_faidx_directory:
    output:
        flag = "Logs/Databases/{db_name_flat}.fa_idx.done"
    input:
        fasta_dir = lambda wc: f"resources/{db_flat_to_path(wc.db_name_flat)}"
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[samtools_faidx_directory]: Indexing all FASTA in {input.fasta_dir}"
    shell:
        r"""
        mkdir -p {input.fasta_dir}
        touch {output.flag}
    
        for fasta in $(find {input.fasta_dir} -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fsa" -o -name "*.fasta" \)); do
            prefix="${{fasta%.*}}"
            echo "Indexing $fasta -> $prefix"
            cmd="samtools faidx $fasta"
            echo "Executing command:\n$cmd\n" >> {output.flag} 
            eval $cmd >> {output.flag} 2>&1
            date -I >> {output.flag}
        done    
        """

##### ALIGNMENT DATABASES #####

rule setup_Ecoli_alignment_db:
    conda:
        config["analysis_settings"]["Escherichia_coli_db"]["yaml"]
    params:
        db_prefix = 'ecoligenes'
    output:
        database = directory(f"{database_path}/{config['analysis_settings']['Escherichia_coli_db']['database']}")
    log:
        stdout = f'Logs/Databases/setup_Ecoli_alignment_db.log'
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
        database = directory(f'{database_path}/{config["analysis_settings"]["Clostridioides_difficile_db"]["database"]}/Toxin')
    params:
        accession_loci = "AM180355.1:tcdA,tcdB,tcdC AF271719.1:cdtA,cdtB",
        db_toxin = "Cdiff_Toxin"
    log:
        stdout = f'Logs/Databases/setup_CdiffToxin.log'
    message:
        "[setup_CdiffToxin]: Setting up C. difficile toxin database"
    shell:
        r"""
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

rule setup_Senterica_alignment_db:
    conda:
        config["analysis_settings"]["Salmonella_enterica_db"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["Salmonella_enterica_db"]["database"]}')
    params:
        MLST_scheme_contigs = ["aroC", "dnaN", "hemD", "hisD", "purE", "sucA", "thrA"],
        MLST_reference = "MLST_Achtman_ref"
    log:
        stdout = f'Logs/Databases/setup_SEnterica.log'
    message:
        "[setup_SEnterica]: Downloading 7 genes for MLST scheme"
    shell:
        """
        set -euo pipefail
        mkdir -p {output.database}/contigs

        for Contig in {params.MLST_scheme_contigs}; do
            echo "Downloading $Contig..." >> {log.stdout}
            cmd="wget https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/${{Contig}}.fasta.gz -O {output.database}/contigs/${{Contig}}.fasta.gz"
            echo -e \"\nExecuting command:\n$cmd\n\" >> {log.stdout}
            eval $cmd >> {log.stdout} 2>&1
            
            cmd="gunzip -d {output.database}/contigs/${{Contig}}.fasta.gz"
            echo -e \"\nExecuting command:\n$cmd\n\" >> {log.stdout}
            eval $cmd >> {log.stdout} 2>&1
        done

        cmd="wget https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/{params.MLST_reference}.fasta -O {output.database}/{params.MLST_reference}.fasta"
        echo -e \"\nExecuting command:\n$cmd\n\" >> {log.stdout}
        eval $cmd >> {log.stdout} 2>&1        

        date -I > {output.database}/creation.date
        """

######### DATABASE REPEAT PATTERNS FOR CDIFF

species_key = config["analysis_settings"]["Clostridioides_difficile_db"]["species_config"]
repeat_list = species_configs[species_key]["analyses_to_run"]["Repeat_identifier"]["repeats"].split()
combo_list = species_configs[species_key]["analyses_to_run"]["Repeat_identifier"]["combos"].split()
TR_repeat_sequences_str = " ".join(repeat_list)
TR_repeat_types_str = " ".join(repeat_list + combo_list)

trst_output_repeat_fa = expand(
    f'{database_path}/{config["analysis_settings"]["Clostridioides_difficile_db"]["database"]}/TRST/{{repeat}}_repeat_sequences.fa',
    repeat=repeat_list
)

rule setup_CdiffTRST:
    conda:
        config["analysis_settings"]["Clostridioides_difficile_db"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["Clostridioides_difficile_db"]["database"]}/TRST'),
        repeat_fa = trst_output_repeat_fa
    params:
        TR_repeat_sequences = TR_repeat_sequences_str,
        TR_repeat_types = TR_repeat_types_str,
    log:
        stdout = f'Logs/Databases/setup_CdiffTRST.log'
    message:
        "[setup_CdiffTRST]: Downloading TRST repeat sequences and types"
    shell:
        r"""
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

rule setup_PlasmidFinder:
    conda:
        config["analysis_settings"]["plasmidfinder"]["yaml"]
    output: 
        database = directory(f'{database_path}/{config["analysis_settings"]["plasmidfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_PlasmidFinder.log'
    message:
        "[setup_PlasmidFinder]: Setting up PlasmidFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git {output.database} > {log.stdout} 2>&1
        
        date -I > {output.database}/creation.date
        """


rule setup_ResFinder:
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["resfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_ResFinder.log'
    message:
        "[setup_ResFinder]: Setting up ResFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git {output.database} > {log.stdout} 2>&1
        
        date -I > {output.database}/creation.date
        """


rule setup_PointFinder:
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["pointfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_PointFinder.log'
    message:
        "[setup_PointFinder]: Setting up PointFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git {output.database} > {log.stdout} 2>&1

        date -I > {output.database}/creation.date
        """

rule setup_DisinFinder:
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["disinfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_DisinFinder.log'
    message:
        "[setup_DisinFinder]: Setting up DisinFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/disinfinder_db.git {output.database} > {log.stdout} 2>&1

        date -I > {output.database}/creation.date
        """


rule setup_VirulenceFinder:
    conda:
        config["analysis_settings"]["virulencefinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["virulencefinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_VirulenceFinder.log'
    message:
        "[setup_VirulenceFinder]: Setting up VirulenceFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git {output.database} > {log.stdout} 2>&1
       
        date -I > {output.database}/creation.date
        """


rule setup_SerotypeFinder:
    conda:
        config["analysis_settings"]["serotypefinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["serotypefinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_SerotypeFinder.log'
    message:
        "[setup_SerotypeFinder]: Setting up SerotypeFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/serotypefinder_db.git {output.database} > {log.stdout} 2>&1

        date -I > {output.database}/creation.date
        """


rule setup_AMRFinder:
    conda:
        config["analysis_settings"]["amrfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["amrfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_AMRFinder.log'
    message:
        "[setup_AMRFinder]: Setting up AMRFinderPlus database"
    shell:
        """
        cmd="amrfinder_update --database $(dirname {output.database}) --force_update"
            
        echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        date -I > {output.database}/creation.date
        """

rule update_MLST:
    conda:
        config["analysis_settings"]["mlst"]["yaml"]
    output:
        datefile = f'{database_path}/mlst/creation.date'
    log:
        stdout = f'Logs/Databases/update_MLST.log'
    message:
        "[update_MLST]: Updating MLST databases."
    shell:
        """
        DIR=$(which mlst)
        MLSTDIR="$DIR/../../db/pubmlst"

        mlst-download_pub_mlst -d $MLSTDIR  > {log.stdout} 2>&1
        mlst-make_blast_db >> {log.stdout} 2>&1 && date -I > {output.datefile}
        """

# Seqsero database is contained within the github and requires github specific installation to utilize it
rule setup_SeqSero2:
    conda:
        config["analysis_settings"]["seqsero2"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["seqsero2"]["database"]}'),
        flag = "Logs/Databases/seqsero2_db.done"
    log:
        stdout = f'Logs/Databases/SeqSero2.log'
    message:
        "[setup_SeqSero2]: Creating SeqSero2 environment and marker directory"
    shell:
        r"""
        set -euo pipefail

        # Clone SeqSero2 only if not already present
        if [ ! -d "SeqSero2" ]; then
            cmd="git clone https://github.com/denglab/SeqSero2.git {output.database}"
            echo -e "Executing command:\n$cmd\n"  >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
        else
            echo "[setup_SeqSero2]: SeqSero2 already cloned, skipping"  >> {log.stdout} 2>&1
        fi

        cmd="python3 -m pip install {output.database}"
        echo -e "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo '[setup_SeqSero2]: Installed SeqSero2 in conda environment' >> {log.stdout} 2>&1

        touch {output.flag}

        echo '[setup_SeqSero2]: Created SeqSero2 completed flag {output.flag}' >> {log.stdout} 2>&1

        date -I > {output.database}/creation.date
        """

rule setup_all_databases:
    input:
        rules.setup_Ecoli_alignment_db.output.database,
        rules.setup_CdiffToxin.output.database,
        rules.setup_Senterica_alignment_db.output.database,
        rules.setup_PlasmidFinder.output.database,
        rules.setup_ResFinder.output.database,
        rules.setup_PointFinder.output.database,
        rules.setup_DisinFinder.output.database,
        rules.setup_VirulenceFinder.output.database,
        rules.setup_AMRFinder.output.database,
        rules.setup_CdiffTRST.output.database,
        rules.setup_SeqSero2.output.database
