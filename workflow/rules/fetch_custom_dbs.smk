rule fetch_genbank:
    params:
        metafile = "%s/{database}_genbank_metafile.tsv" %metadata_path,
        merge = 500
    output:
        fasta = "%s/custom/{database,[^/]+}.fasta" % database_path,
        bed = "%s/custom/{database,[^/]+}.bed6" % database_path,
    conda:
        "../envs/fetch.yaml"
    log:
        stdout = 'Logs/Databases/fetch_genbank_{database}.log'
    message:
        "[fetch_genbank]: Fetching {wildcards.database} from Genbank"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.fasta})

        cmd="python workflow/scripts/genbank_fetcher.py --metafile {params.metafile} --bed {output.bed} --fasta {output.fasta} --merge {params.merge} --append"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """


rule fetch_type_repeat_sequence:
    output:
        seq = "%s/custom/type_repeats/{TR}.fasta" %database_path
    conda:
        "../envs/fetch.yaml"
    log:
        stdout = "Logs/Databases/fetch_type_repeat_sequences_{TR}.log"
    message:
        "[fetch_type_repeat_sequences]: Downloading Type Repeat Sequence Type sequences"
    shell:
        """
        mkdir -p $(dirname {output.seq})

        cmd="curl -fSL https://raw.githubusercontent.com/ssi-dk/cdiff_fbi/refs/heads/raah_dev/db/TRST/{wildcards.TR}_repeat_sequences.fa -o {output.seq}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """


rule fetch_type_repeat_metadata:
    output:
        meta = "%s/custom/type_repeats/{TR}.txt" %database_path
    conda:
        "../envs/fetch.yaml"
    log:
        stdout = "Logs/Databases/fetch_type_repeat_metadata_{TR}.log"
    message:
        "[fetch_type_repeat_metadata]: Downloading Type Repeat Sequence Type metadata"
    shell:
        """
        mkdir -p $(dirname {output.meta})

        cmd="curl -fSL https://raw.githubusercontent.com/ssi-dk/cdiff_fbi/refs/heads/raah_dev/db/TRST/{wildcards.TR}_repeat_types.txt -o {output.meta}"
        
        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """


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


rule fetch_Senterica_Scheme:
  output:
    source = "%s/custom/SalmonellaAchtman7GeneMLST.fasta" %database_path,
    profile = "%s/custom/SalmonellaAchtman7GeneMLST.txt" %database_path
  conda:
    "../envs/fetch.yaml"
  log:
    stdout = "Logs/Databases/setup_senterica.log"
  message:
    "[fetch_Senterica_Scheme]: Downloading Achtman 7 Gene MLST scheme for Salmonella Enterica"
  shell:
    """
    mkdir -p $(dirname {output.source})
    cmd="curl https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/MLST_Achtman_ref.fasta -o {output.source}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1

    cmd="curl https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/profiles.list.gz | gunzip -c > {output.profile}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


rule fetch_Senterica_Serovar:
  output:
    source = "%s/custom/Senterica_serovar.txt" %database_path
  conda:
    "../envs/fetch.yaml"
  log:
    stdout = "Logs/Databases/setup_senterica.log"
  message:
    "[fetch_Senterica_Scheme]: Downloading Achtman 7 Gene MLST scheme for Salmonella Enterica"
  shell:
    """
    mkdir -p $(dirname {output.source})
    cmd="curl -fSL https://raw.githubusercontent.com/phac-nml/sistr_cmd/master/sistr/data/serovar-list.txt -o {output.source}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


rule fetch_chtyper_db:
    output:
        source = "%s/custom/fumCH_db.fasta" %database_path
    conda:
        "../envs/fetch.yaml"
    log:
        stdout = "Logs/Databases/setup_chtyper_database.log"
    message:
        "[fetch_chtyper_db]: Downloading custom database for CHtyper"
    shell:
        """
        outdir=$(dirname {output.source})
        mkdir -p $outdir
    
        cmd="curl https://bitbucket.org/genomicepidemiology/chtyper_db/raw/654ca48d250e0a69c6c06b4be5a96d807b23f806/fimH.fsa -o $outdir/fimH.fsa"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    
        cmd="curl https://bitbucket.org/genomicepidemiology/chtyper_db/raw/654ca48d250e0a69c6c06b4be5a96d807b23f806/fumC.fsa -o $outdir/fumC.fsa"
        eval $cmd >> {log.stdout} 2>&1 

        cmd="cat $outdir/fimH.fsa $outdir/fumC.fsa > {output.source}"
        eval $cmd >> {log.stdout} 2>&1
        """

# Place holder rule until we have an online repo for all dbs
# We store momentarily in Dataset/databases
rule fetch_blast_database:
    conda:
        "../envs/blast.yaml"
    input: 
        source = "%s/{database}.fasta" %temp_storage_path
    output:
        source = "%s/custom/blast/{database}.fasta" %database_path
    log:
        stdout = 'Logs/Databases/setup_{database}.log'
    message:
        "[Fetch {wildcards.database} Blast database]: Setting up the {wildcards.database} database from the temporary storage folder"
    shell:
        """
        cmd="cp {input.source} {output.source}"
            
        echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """
