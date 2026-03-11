rule fetch_genbank:
    params:
        metafile = "%s/{database}_genbank_metafile.tsv" %TARGET_SCREENING_DIR,
        merge = 500
    output:
        fasta = "%s/custom/{database,[^/]+}.fasta" % database_dir,
        bed = "%s/custom/{database,[^/]+}.bed6" % database_dir,
        version_db = "%s/custom/{database,[^/]+}_version.txt" %database_dir
    conda:
        ENVS_DIR / "fetch.yaml"
    log:
        stdout = "%s/Databases/fetch_genbank_{database}.log" %logdir
    message:
        "[fetch_genbank]: Fetching {wildcards.database} from Genbank"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.fasta})

        # 1) Run the genbank fetcher
        cmd="python {SCRIPTS_DIR}/genbank_fetcher.py --metafile {params.metafile} --bed {output.bed} --fasta {output.fasta} --merge {params.merge} --append"

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
        seq = "%s/custom/type_repeats/{TR}.fasta" %database_dir,
        version_db = "%s/custom/type_repeats/{TR}_version.txt" % database_dir
    conda:
        ENVS_DIR / "fetch.yaml"
    log:
        stdout = "%s/Databases/fetch_type_repeat_sequences_{TR}.log" %logdir
    message:
        "[fetch_type_repeat_sequences]: Downloading Type Repeat Sequence Type sequences"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.seq})

        rel_path="clostridioides_difficile/type_repeats/{wildcards.TR}.fasta"
        rel_ver="clostridioides_difficile/type_repeats/{wildcards.TR}_version.txt"

        fasta_url="https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/main/$rel_path"
        ver_url="https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/main/$rel_ver"

        cmd_fasta="curl -fSL $fasta_url -o {output.seq}"
        cmd_ver="curl -fSL $ver_url -o {output.version_db}"

        echo "Executing command:\n$cmd_fasta\n$cmd_ver\n" > {log.stdout}
        eval "$cmd_fasta" >> {log.stdout} 2>&1
        eval "$cmd_ver"   >> {log.stdout} 2>&1
        """


rule fetch_type_repeat_metadata:
    output:
        meta = "%s/custom/type_repeats/{TR}.txt" %database_dir
    conda:
        ENVS_DIR / "fetch.yaml"
    log:
        stdout = "%s/Databases/fetch_type_repeat_metadata_{TR}.log" %logdir
    message:
        "[fetch_type_repeat_metadata]: Downloading Type Repeat Sequence Type metadata"
    shell:
        """
        mkdir -p $(dirname {output.meta})

        cmd="curl -fSL https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/refs/heads/main/clostridioides_difficile/type_repeats/{wildcards.TR}.txt -o {output.meta}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """

rule fetch_ecoligenes:
    output:
        source = "%s/custom/ecoligenes.fasta" %database_dir,
        version_db = "%s/custom/ecoligenes_version.txt" % database_dir
    conda:
        ENVS_DIR / "fetch.yaml"
    log:
        stdout = "%s/Databases/setup_ecoligenes_ecoligenes.log" %logdir
    message:
        "[fetch_ecoligenes]: Downloading custom database ecoligenes"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.source})
        
        rel_path="escherichia_coli/ecoligenes.fasta"
        rel_ver="escherichia_coli/ecoligenes_version.txt"
        
        fasta_url="https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/main/$rel_path"
        ver_url="https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/main/$rel_ver"

        cmd_fasta="curl -fSL $fasta_url -o {output.source}"
        cmd_ver="curl -fSL $ver_url | awk '1' > {output.version_db}"

        echo "Executing command:\n$cmd_fasta\n$cmd_ver\n" > {log.stdout}
        eval "$cmd_fasta" >> {log.stdout} 2>&1
        eval "$cmd_ver"   >> {log.stdout} 2>&1
        """

rule fetch_Senterica_Scheme:
    output:
        source = "%s/custom/SalmonellaAchtman7GeneMLST.fasta" %database_dir,
        profile = "%s/custom/SalmonellaAchtman7GeneMLST.txt" %database_dir,
        version_db = "%s/custom/SalmonellaAchtman7GeneMLST_version.txt" % database_dir
    conda:
        ENVS_DIR / "fetch.yaml"
    log:
        stdout = "%s/Databases/setup_senterica_mlst_scheme.log" %logdir
    message:
        "[fetch_Senterica_Scheme]: Downloading Achtman 7 Gene MLST scheme for Salmonella Enterica"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.source})
                
        fasta_url="https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/MLST_Achtman_ref.fasta"
        profile_url="https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/profiles.list.gz"

        fasta_cmd="curl -fSL $fasta_url -o {output.source}"
        profile_cmd="curl -fSL $profile_url -o {output.profile}"

        echo "Executing command:\n$fasta_cmd\n" > {log.stdout} 2>&1
        eval $fasta_cmd >> {log.stdout} 2>&1

        echo "Executing command:\n$profile_cmd\n" > {log.stdout} 2>&1
        eval $profile_cmd >> {log.stdout} 2>&1
        
        # 2) Get ETag (as a clean value) and make version file
        etag_cmd="curl -sI $fasta_url | sed -n 's/^etag: //Ip' | tr -d '\\r' | tr -d '\\042'"
        date_cmd="date -I"

        echo -e "Executing command:\n$etag_cmd\n$date_cmd\n" >> {log.stdout}

        etag_str="$(eval "$etag_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        # Fallback if no ETag is present for some reason
        if [ -z "$etag_str" ]; then
            etag_str="no_etag"
        fi

        # Build version ID. If you DON'T want the ETag at all, set version_str="SalmonellaAchtman7GeneMLST"
        version_str="SalmonellaAchtman7GeneMLST_$etag_str"

        # Write "<id>\t<download_date>" to the version file
        printf '%s\t%s\n' "$version_str" "$date_str" > {output.version_db}
        """

rule fetch_Senterica_Serovar:
    output:
        source = "%s/custom/Senterica_serovar.txt" % database_dir,
        version_db = "%s/custom/Senterica_serovar_version.txt" % database_dir
    conda:
        ENVS_DIR / "fetch.yaml"
    log:
        stdout = "%s/Databases/setup_senterica_sistr.log" %logdir
    message:
        "[fetch_Senterica_Serovar]: Downloading SISTR serovar list"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.source})

        list_url="https://raw.githubusercontent.com/phac-nml/sistr_cmd/master/sistr/data/serovar-list.txt"

        # 1) Download the serovar list
        cmd_fasta="curl -fSL $list_url -o {output.source}"

        echo "Executing command:\n$cmd_fasta\n" > {log.stdout}
        eval "$cmd_fasta" >> {log.stdout} 2>&1

        # 2) Get ETag (as a clean value) and make version file
        # sed -n 's/^etag: //Ip' greps case insensitive etag and replace with nothing
        # tr -d \\r deletes the characters since some HTTP headers end lines with \r\n
        # tr -d 042 deletes the octal character for double quotes
        etag_cmd="curl -sI $list_url | sed -n 's/^etag: //Ip' | tr -d '\\r' | tr -d '\\042'"
        date_cmd="date -I"

        echo -e "Executing command:\n$etag_cmd\n$date_cmd\n" >> {log.stdout}

        etag_str="$(eval "$etag_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        # Fallback if no ETag is present for some reason
        if [ -z "$etag_str" ]; then
            etag_str="no_etag"
        fi

        # Build version ID. If you DON'T want the ETag at all, set version_str="sistr_serovar_list"
        version_str="sistr_serovar_list_$etag_str"

        # Write "<id>\t<download_date>" to the version file
        printf '%s\t%s\n' "$version_str" "$date_str" > {output.version_db}
        """

rule setup_LREfinder:
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    params:
        prefix = "%s/custom/" %database_dir,
        dbdir = "%s/custom/elmDB/" %database_dir,
    output:
        source = "%s/custom/elmDB.fasta" %database_dir,
        version_db = "%s/custom/elmDB_version.txt" % database_dir
    log:
        stdout = "%s/Databases/LREfinder_db.log" %logdir
    message:
        "[setup_LREfinder]: Setting up LREfinder database"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.source})

        sequence_url="https://bitbucket.org/genomicepidemiology/lre-finder/raw/fac445d190853cc90c1aed392a55102fe9df4376/elmDB.tar.gz"

        # 1) download raw sequence
        curl -fSL $sequence_url --output - | tar -xzvf - -C {params.prefix}
        mv {params.dbdir}elm.fsa {output.source}
        rm -rf {params.dbdir}

        # 2) download version from etag
        etag_cmd="curl -sI $sequence_url | sed -n 's/^etag: //Ip' | tr -d '\\r' | tr -d '\\042'"
        date_cmd="date -I"

        echo -e "Executing command:\n$etag_cmd\n$date_cmd\n" >> {log.stdout}

        etag_str="$(eval "$etag_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        # Fallback if no ETag is present for some reason
        if [ -z "$etag_str" ]; then
            etag_str="no_etag"
        fi
        
        # Build version ID. If you DON'T want the ETag at all, set version_str="sistr_serovar_list"
        version_str="LREfinder_elmDB_$etag_str"

        # Write "<id>\t<download_date>" to the version file
        printf '%s\t%s\n' "$version_str" "$date_str" > {output.version_db}
        """

rule fetch_chtyper_db:
    output:
        source = "%s/custom/fumCH_db.fasta" %database_dir,
        version_db = "%s/custom/fumCH_db_version.txt" % database_dir
    conda:
        ENVS_DIR / "fetch.yaml"
    log:
        stdout = "%s/Databases/setup_chtyper_database.log" %logdir
    message:
        "[fetch_chtyper_db]: Downloading custom database for CHtyper"
    shell:
        """
        set -euo pipefail
        outdir=$(dirname {output.source})
        mkdir -p $outdir
        
        fimH_url="https://bitbucket.org/genomicepidemiology/chtyper_db/raw/654ca48d250e0a69c6c06b4be5a96d807b23f806/fimH.fsa"
        fumC_url="https://bitbucket.org/genomicepidemiology/chtyper_db/raw/654ca48d250e0a69c6c06b4be5a96d807b23f806/fumC.fsa"

        # 1) Download the serovar list
        cmd_fimH="curl -fSL $fimH_url -o $outdir/fimH.fsa"
        cmd_fumC="curl -fSL $fumC_url -o $outdir/fumC.fsa"

        echo "Executing command:\n$cmd_fimH\n" > {log.stdout}
        eval "$cmd_fimH" >> {log.stdout} 2>&1
        echo "Executing command:\n$cmd_fumC\n" >> {log.stdout}
        eval "$cmd_fumC" >> {log.stdout} 2>&1

        # 2) Get the etag versions
        etag_cmd="curl -sI $fimH_url | sed -n 's/^etag: //Ip' | tr -d '\\r' | tr -d '\\042'"
        date_cmd="date -I"

        echo -e "Executing command:\n$etag_cmd\n$date_cmd\n" >> {log.stdout}

        etag_str1="$(eval "$etag_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        # 3) create final database
        cmd="cat $outdir/fimH.fsa $outdir/fumC.fsa > {output.source}"
        eval $cmd >> {log.stdout} 2>&1

        # Write "<accessions>\t<date>" to the version file
        printf '%s%s\t%s\n' "chtyper_" "$etag_str1" "$date_str" > {output.version_db}
        """

# Place holder rule until we have an online repo for all dbs
# We store momentarily in Dataset/databases
rule fetch_custom_blast_database:
    conda:
        ENVS_DIR / "blast.yaml"
    output:
        source = "%s/custom/blast/OXAndm.fasta" %database_dir,
        version_db = "%s/custom/blast/OXAndm_version.txt" % database_dir
    log:
        stdout = "%s/Databases/setup_OXAndm.log" %logdir
    message:
        "[fetch_custom_blast_database]: Downloading custom database OXAndm"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.source})
        
        rel_path="escherichia_coli/OXAndm.fasta"
        rel_ver="escherichia_coli/OXAndm_version.txt"
        
        fasta_url="https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/main/$rel_path"
        ver_url="https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/main/$rel_ver"

        cmd_fasta="curl -fSL $fasta_url -o {output.source}"
        cmd_ver="curl -fSL $ver_url  | awk '1' > {output.version_db}"

        echo "Executing command:\n$cmd_fasta\n$cmd_ver\n" > {log.stdout}
        eval "$cmd_fasta" >> {log.stdout} 2>&1
        eval "$cmd_ver"   >> {log.stdout} 2>&1
        """
