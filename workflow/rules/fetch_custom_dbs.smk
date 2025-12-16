rule fetch_genbank:
    params:
        metafile = "%s/{database}_genbank_metafile.tsv" %metadata_path,
        merge = 500
    output:
        fasta = "%s/custom/{database,[^/]+}.fasta" % database_path,
        bed = "%s/custom/{database,[^/]+}.bed6" % database_path,
        version_db = "%s/custom/{database,[^/]+}_version.txt" %database_path
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

        # 1) Run the genbank fetcher
        cmd="python workflow/scripts/genbank_fetcher.py --metafile {params.metafile} --bed {output.bed} --fasta {output.fasta} --merge {params.merge} --append"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1

        # 2) Make version file with one '_'-separated line of unique accessions and starting with genbank_db_
        
        version_cmd="tail -n +2 {params.metafile} | cut -f1 | sort -u | paste -sd '_' - | sed 's/^/genbank_db_/'"
        date_cmd="date -I"

        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        # Write "<accessions>\t<date>" to the version file
        printf '%s\t%s\n' "$version_str" "$date_str" > {output.version_db}
        """


rule fetch_type_repeat_sequence:
    output:
        seq = "%s/custom/type_repeats/{TR}.fasta" %database_path,
        version_db = "%s/custom/type_repeats/{TR}_version.txt" % database_path
    conda:
        "../envs/fetch.yaml"
    log:
        stdout = "Logs/Databases/fetch_type_repeat_sequences_{TR}.log"
    message:
        "[fetch_type_repeat_sequences]: Downloading Type Repeat Sequence Type sequences"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.seq})

        rel_path="clostridioides_difficile/type_repeats/{wildcards.TR}.fasta"
        rel_ver="clostridioides_difficile/type_repeats/{wildcards.TR}_version.txt"

        fasta_url="https://raw.githubusercontent.com/RAHenriksen/ssi_analysis_utility_db/main/$rel_path"
        ver_url="https://raw.githubusercontent.com/RAHenriksen/ssi_analysis_utility_db/main/$rel_ver"

        cmd_fasta="curl -fSL $fasta_url -o {output.seq}"
        cmd_ver="curl -fSL $ver_url -o {output.version_db}"

        echo "Executing command:\n$cmd_fasta\n$cmd_ver\n" > {log.stdout}
        eval "$cmd_fasta" >> {log.stdout} 2>&1
        eval "$cmd_ver"   >> {log.stdout} 2>&1
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

        cmd="curl -fSL https://raw.githubusercontent.com/RAHenriksen/ssi_analysis_utility_db/refs/heads/main/clostridioides_difficile/type_repeats/{wildcards.TR}.txt -o {output.meta}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """


rule fetch_ecoligenes:
    output:
        source = "%s/custom/ecoligenes.fasta" %database_path,
        version_db = "%s/custom/ecoligenes_version.txt" % database_path
    conda:
        "../envs/fetch.yaml"
    log:
        stdout = "Logs/Databases/setup_ecoligenes_ecoligenes.log"
    message:
        "[fetch_ecoligenes]: Downloading custom database ecoligenes"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.source})
        
        rel_path="escherichia_coli/ecoligenes.fasta"
        rel_ver="escherichia_coli/ecoligenes_version.txt"
        
        fasta_url="https://raw.githubusercontent.com/RAHenriksen/ssi_analysis_utility_db/main/$rel_path"
        ver_url="https://raw.githubusercontent.com/RAHenriksen/ssi_analysis_utility_db/main/$rel_ver"

        cmd_fasta="curl -fSL $fasta_url -o {output.source}"
        cmd_ver="curl -fSL $ver_url -o {output.version_db}"

        echo "Executing command:\n$cmd_fasta\n$cmd_ver\n" > {log.stdout}
        eval "$cmd_fasta" >> {log.stdout} 2>&1
        eval "$cmd_ver"   >> {log.stdout} 2>&1
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

rule setup_LREfinder:
    conda:
        "../envs/kmeraligner.yaml"
    output:
        database = "%s/custom/elmDB.fasta" %database_path
    log:
        stdout = 'Logs/Databases/LREfinder_db.log'
    message:
        "[setup_LREfinder]: Setting up LREfinder database"
    shell:
        """
        outdir=$(dirname {output.database})
        mkdir -p $outdir
        curl -fSL https://bitbucket.org/genomicepidemiology/lre-finder/raw/fac445d190853cc90c1aed392a55102fe9df4376/elmDB.tar.gz --output elmDB.tar.gz 
        tar -xvf elmDB.tar.gz
        mv elmDB/elm.fsa {output.database}
        rm -r elmDB/ elmDB.tar.gz
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
