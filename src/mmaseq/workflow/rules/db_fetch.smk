rule fetch_genbank:
    params:
        metafile = f"{SCREENING_DIR}/{{database}}_genbank_metafile.tsv",
        merge = 500
    output:
        fasta = f"{database_dir}/custom/{{database}}.fasta",
        bed = f"{database_dir}/custom/{{database}}.bed6",
        version_db = f"{database_dir}/custom/{{database}}_version.txt"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/Databases/fetch_genbank_{{database}}.log"
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
        seq = f"{database_dir}/type_repeats/{{TR}}.fasta",
        version_db = f"{database_dir}/type_repeats/{{TR}}_version.txt"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/Databases/fetch_type_repeat_sequence_{{TR}}.log"
    message:
        "[fetch_type_repeat_sequence]: Downloading Type Repeat Sequence Type sequences"
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
        meta = f"{database_dir}/type_repeats/{{TR}}.txt"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/Databases/fetch_type_repeat_metadata_{{TR}}.log"
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
        source = f"{database_dir}/custom/ecoligenes.fasta",
        version_db = f"{database_dir}/custom/ecoligenes_version.txt"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/Databases/fetch_ecoligenes.log"
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

rule fetch_senterica_scheme:
    output:
        source = f"{database_dir}/custom/SalmonellaAchtman7GeneMLST.fasta",
        profile = f"{database_dir}/custom/SalmonellaAchtman7GeneMLST.txt",
        version_db = f"{database_dir}/custom/SalmonellaAchtman7GeneMLST_version.txt"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/Databases/fetch_senterica_scheme.log"
    message:
        "[fetch_senterica_scheme]: Downloading Achtman 7 Gene MLST scheme for Salmonella Enterica"
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

rule fetch_senterica_serovar:
    output:
        source = f"{database_dir}/custom/Senterica_serovar.txt",
        version_db = f"{database_dir}/custom/Senterica_serovar_version.txt"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/Databases/fetch_senterica_serovar.log"
    message:
        "[fetch_senterica_serovar]: Downloading SISTR serovar list"
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


rule fetch_chtyper_db:
    output:
        source = f"{database_dir}/custom/fumCH_db.fasta",
        version_db = f"{database_dir}/custom/fumCH_db_version.txt"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/Databases/fetch_chtyper_db.log"
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
    output:
        source = f"{database_dir}/custom/blast/OXAndm.fasta",
        version_db = f"{database_dir}/custom/blast/OXAndm_version.txt"
    conda:
        ENVS_DIR / "blast.yaml"
    log:
        stdout = f"{logdir}/Databases/fetch_custom_blast_database.log"
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


# Place holder rule until we have an online repo for all dbs
# We store momentarily in Dataset/databases
rule fetch_vancomycin:
    output:
        source = f"{database_dir}/custom/vancomycin.fasta"
    log:
        stdout = f"{logdir}/Databases/fetch_vancomycin.log"
    message:
        "[fetch_vancomycin]: Downloading custom vancomycin database"
    shell:
        """
        set -euo pipefail
        OUTDIR=$(dirname {output.source})
        mkdir -p $OUTDIR
                
        fasta_url="https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/refs/heads/main/resistance/vancomycin.fasta"

        cmd="curl -fSL $fasta_url -o {output.source}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval "$cmd" >> {log.stdout} 2>&1
        """


rule fetch_vancomycin_operon:
    output:
        source = f"{database_dir}/custom/vancomycinOperon.fasta"
    log:
        stdout = f"{logdir}/Databases/fetch_vancomycin_operon.log"
    message:
        "[fetch_vancomycin_operon]: Downloading custom vancomycin operon database"
    shell:
        """
        set -euo pipefail
        OUTDIR=$(dirname {output.source})
        mkdir -p $OUTDIR
                
        fasta_url="https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/refs/heads/main/resistance/vancomycinOperon.fasta"

        cmd="curl -fSL $fasta_url -o {output.source}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval "$cmd" >> {log.stdout} 2>&1
        """