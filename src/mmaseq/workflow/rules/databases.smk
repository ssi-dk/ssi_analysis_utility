### Processing ###

rule setup_PlasmidFinder:
    output: 
        database = directory("%s/plasmidfinder_db" %database_dir),
        version_db = "%s/plasmidfinder_db/PlasmidFinder_version.txt" %database_dir
    conda:
        ENVS_DIR / "plasmidfinder.yaml"
    log:
        stdout = "%s/Databases/setup_PlasmidFinder.log" %logdir
    message:
        "[setup_PlasmidFinder]: Setting up PlasmidFinder database"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.database})

        # 1) Download database
        git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[plasmidfinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the plasmidfinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        # 2) create version file with date
        version_cmd="cat \"{output.database}/VERSION\""
        date_cmd="date -I"
        
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s\t%s\n' "PlasmidFinder_" "$version_str" "$date_str" > {output.version_db}
        """

rule setup_ResFinder:
    output:
        database = directory("%s/resfinder_db" %database_dir),
        version_db = "%s/resfinder_db/ResFinder_version.txt" %database_dir
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = "%s/Databases/setup_ResFinder.log" %logdir
    message:
        "[setup_ResFinder]: Setting up ResFinder database"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.database})

        # 1) Download database
        git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git {output.database} > {log.stdout} 2>&1
        
        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"

            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[resfinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the resfinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        # 2) create version file with date
        version_cmd="cat \"{output.database}/VERSION\""
        date_cmd="date -I"
        
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s\t%s\n' "ResFinder_" "$version_str" "$date_str" > {output.version_db}
        """

rule setup_PointFinder:
    output:
        database = directory("%s/pointfinder_db" %database_dir),
        version_db = "%s/pointfinder_db/PointFinder_version.txt" %database_dir
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = "%s/Databases/setup_PointFinder.log" %logdir
    message:
        "[setup_PointFinder]: Setting up PointFinder database"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.database})

        # 1) Download database
        git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -type f -name '*.fsa'); do
            idx_prefix="${{fasta%.*}}"

            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1

            if [ -z $idx_prefix.comb.b ]; then
                echo '[pointfinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the pointfinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        # 2) create version file with date
        version_cmd="cat \"{output.database}/VERSION\""
        date_cmd="date -I"
        
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s\t%s\n' "PointFinder_" "$version_str" "$date_str" > {output.version_db}
        """

rule setup_DisinFinder:
    output:
        database = directory("%s/disinfinder_db" %database_dir),
        version_db = "%s/disinfinder_db/DisinFinder_version.txt" %database_dir
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = "%s/Databases/setup_DisinFinder.log" %logdir
    message:
        "[setup_DisinFinder]: Setting up DisinFinder database"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.database})

        # 1) Download database
        git clone https://bitbucket.org/genomicepidemiology/disinfinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[disinfinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the disinfinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        # 2) create version file with date
        version_cmd="cat \"{output.database}/VERSION\""
        date_cmd="date -I"
        
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s\t%s\n' "DisinFinder_" "$version_str" "$date_str" > {output.version_db}
        """

rule setup_VirulenceFinder:
    output:
        database = directory("%s/virulencefinder_db" %database_dir),
        version_db = "%s/virulencefinder_db/VirulenceFinder_version.txt" %database_dir
    conda:
        ENVS_DIR / "virulencefinder.yaml"
    log:
        stdout = "%s/Databases/setup_VirulenceFinder.log" %logdir
    message:
        "[setup_VirulenceFinder]: Setting up VirulenceFinder database"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.database})

        # 1) Download database
        git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git {output.database} > {log.stdout} 2>&1
       
        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[virulencefinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the virulencefinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        # 2) create version file with date
        version_cmd="cat \"{output.database}/VERSION\""
        date_cmd="date -I"
        
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s\t%s\n' "VirulenceFinder_" "$version_str" "$date_str" > {output.version_db}
        """

rule setup_SerotypeFinder:
    output:
        database = directory("%s/serotypefinder_db" %database_dir),
        version_db = "%s/serotypefinder_db/SerotypeFinder_version.txt" %database_dir
    conda:
        ENVS_DIR / "serotypefinder.yaml"
    log:
        stdout = "%s/Databases/setup_SerotypeFinder.log" %logdir
    message:
        "[setup_SerotypeFinder]: Setting up SerotypeFinder database"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.database})

        # 1) Download database
        git clone https://bitbucket.org/genomicepidemiology/serotypefinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[serotypefinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the serotypefinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        # 2) create version file with date
        version_cmd="cat \"{output.database}/VERSION\""
        date_cmd="date -I"
        
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s\t%s\n' "SerotypeFinder_" "$version_str" "$date_str" > {output.version_db}
        """

rule setup_Spatyper:
    output:
        database = directory("%s/spatyper_db" %database_dir)
    log:
        stdout = "%s/Databases/setup_Spatyper.log" %logdir
    message:
        "[Setup Spatyper]: Setting up SerotypeFinder database"
    shell:
        """
        cmd="git clone https://bitbucket.org/genomicepidemiology/spatyper_db.git {output.database}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        date -I > {output.database}/creation.date
        """


rule setup_AMRFinder:
    output:
        database = directory("%s/amrfinderplus/latest" %database_dir),
        version_db = "%s/amrfinderplus/latest/AMRFinder_version.txt" %database_dir
    conda:
        ENVS_DIR / "amrfinder.yaml"
    log:
        stdout = "%s/Databases/setup_AMRFinder.log" %logdir
    message:
        "[setup_AMRFinder]: Setting up AMRFinderPlus database"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.database})

        # 1) Download database
        cmd="amrfinder_update --database $(dirname {output.database}) --force_update"
            
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # 2) create version file with date
        modified_cmd="cat \"{output.database}/version.txt\""
        version_cmd="amrfinder_update --version"
        date_cmd="date -I"
        
        echo -e "Executing command:\n$modified_cmd\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        modified_str="$(eval "$modified_cmd" 2>> {log.stdout})"
        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s%s%s\t%s\n' "AMRFinder_" "$version_str" "_" "$modified_str" "$date_str" > {output.version_db}
        """


rule setup_kleborate_amrfinder:
    input:
        database = rules.setup_AMRFinder.output.database,
        version_db = rules.setup_AMRFinder.output.version_db
    output:
        version_db = "%s/kleborate/kleborate_version.txt" %database_dir
    conda:
        ENVS_DIR / "kleborate.yaml"
    log:
        stdout = "%s/Databases/setup_Kleborate_AMRFinder.log" %logdir
    message:
        "[Setup_Kleborate_AMRFinder]: Clonig AMRFinder database to Kleborate environment"
    shell:
        """
            BIN=$(which kleborate)
            BINDIR=$(dirname $BIN)
            KLEBDIR="$(dirname $BINDIR)/share/amrfinderplus/data"

            mkdir -p  $KLEBDIR
            database=$(realpath -s {input.database})
            cmd="ln -sf $database $KLEBDIR/"

            echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1

            mkdir -p $(dirname {output.version_db})
            version_database=$(realpath -s {input.version_db})
            cmd="ln -sf $version_database {output.version_db}"
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
        """


rule setup_LREfinder:
    params:
        prefix = "%s/custom/" %database_dir,
        dbdir = "%s/custom/elmDB/" %database_dir,
    output:
        source = "%s/custom/elmDB.fasta" %database_dir,
        version_db = "%s/custom/elmDB_version.txt" % database_dir
    conda:
        ENVS_DIR / "kmeraligner.yaml"
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
        cmd="curl -fSL $sequence_url --output - | tar -xzvf - -C {params.prefix} && mv {params.dbdir}elm.fsa {output.source} && rm -rf {params.dbdir}"

        echo -e "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1


        # 2) download version from etag
        etag_cmd="curl -sI $sequence_url | sed -n 's/^etag: //Ip' | tr -d '\\r' | tr -d '\\042'"
        date_cmd="date -I"

        echo -e "Executing command:\n$etag_cmd\n$date_cmd\n" >> {log.stdout} 2>&1

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


rule setup_custom_kmeraligner_index:
    input:
        source = "%s/custom/{database}.fasta" %database_dir
    params:
        prefix = "%s/kmeraligner/{database}" %database_dir
    output:
        combined_size = "%s/kmeraligner/{database}.comp.b" %database_dir,
        lengths = "%s/kmeraligner/{database}.length.b" %database_dir,
        names = "%s/kmeraligner/{database}.name" %database_dir,
        seqs = "%s/kmeraligner/{database}.seq.b" %database_dir,
        version_db = "%s/kmeraligner/{database}_kmaindex_version.txt" %database_dir
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = "%s/Databases/setup_custom_kmeraligner_index_{database}.log" %logdir
    message:
        "[setup_custom_kmeraligner_index]: Setting up {wildcards.database} database with kmeraligner"
    shell:
        """
            mkdir -p $(dirname {params.prefix})

            cmd="kma index -i {input.source} -o {params.prefix}"

            echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1

            # 2) create version file with date
            version_cmd="kma index -v"
            date_cmd="date -I"
                
            echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

            version_str="$(eval "$version_cmd" 2>> {log.stdout})"
            date_str="$(eval "$date_cmd" 2>> {log.stdout})"

            printf '%s\t%s\n' "$version_str" "$date_str" > {output.version_db}
        """

rule setup_custom_bowtie2_index:
    input:
        source = "%s/custom/{database}.fasta" %database_dir
    params:
        prefix = "%s/bowtie2/{database}" %database_dir
    output:
        bt2_1 = "%s/bowtie2/{database}.1.bt2" %database_dir,
        bt2_2 = "%s/bowtie2/{database}.2.bt2" %database_dir,
        bt2_3 = "%s/bowtie2/{database}.3.bt2" %database_dir,
        bt2_4 = "%s/bowtie2/{database}.4.bt2" %database_dir,
        bt2_1_rev = "%s/bowtie2/{database}.rev.1.bt2" %database_dir,
        bt2_2_rev = "%s/bowtie2/{database}.rev.2.bt2" %database_dir,
        version_db = "%s/bowtie2/{database}_bowtie2index_version.txt" %database_dir
    conda:
        ENVS_DIR / "bowtie2.yaml"
    log:
        stdout = "%s/Databases/setup_custom_bowtie2index_{database}.log" %logdir
    message:
        "[setup_custom_bowtie2_index]: Setting up {wildcards.database} database with bowtie2"
    shell:
        """
            mkdir -p $(dirname {params.prefix})
            cmd="bowtie2-build {input.source} {params.prefix}"

            echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1

            # 2) create version file with date
            version_cmd="bowtie2-build --version | head -n1 | grep -oE '[0-9]+([.][0-9]+)+'"
            date_cmd="date -I"
                
            echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

            version_str="$(eval "$version_cmd" 2>> {log.stdout})"
            date_str="$(eval "$date_cmd" 2>> {log.stdout})"

            printf '%s\t%s\n' "$version_str" "$date_str" > {output.version_db}
        """


rule setup_custom_samtool_index:
    input:
        source = "%s/custom/{database}.fasta" %database_dir
    output:
        source = "%s/samtools/{database}.fasta" %database_dir, 
        index = "%s/samtools/{database}.fasta.fai" %database_dir
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = "%s/Databases/setup_custom_samtool_index_{database}.log" %logdir
    message:
        "[setup_custom_samtool_index]: Setting up {wildcards.database} database with samtools"
    shell:
        """
            outdir=$(dirname {output.source})
            mkdir -p $outdir

            cmd="cp {input.source} {output.source}"

            echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1

            cmd="samtools faidx {output.source} -o {output.index}"

            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1

            date -I > $outdir/creation.date
        """


### Pure downloaders ###

rule fetch_genbank:
    params:
        metafile = "%s/{database}_genbank_metafile.tsv" %SCREENING_DIR,
        merge = 500
    output:
        fasta = "%s/custom/{database,[^/]+}.fasta" % database_dir,
        bed = "%s/custom/{database,[^/]+}.bed6" % database_dir,
        version_db = "%s/custom/{database,[^/]+}_version.txt" %database_dir
    conda:
        ENVS_DIR / "py_utls.yaml"
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
        ENVS_DIR / "py_utls.yaml"
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
        ENVS_DIR / "py_utls.yaml"
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
        ENVS_DIR / "py_utls.yaml"
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
        ENVS_DIR / "py_utls.yaml"
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
        ENVS_DIR / "py_utls.yaml"
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


rule fetch_chtyper_db:
    output:
        source = "%s/custom/fumCH_db.fasta" %database_dir,
        version_db = "%s/custom/fumCH_db_version.txt" % database_dir
    conda:
        ENVS_DIR / "py_utls.yaml"
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


rule fetch_custom_blast_database:
    output:
        source = "%s/custom/blast/OXAndm.fasta" %database_dir,
        version_db = "%s/custom/blast/OXAndm_version.txt" % database_dir
    conda:
        ENVS_DIR / "blast.yaml"
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


rule fetch_vancomycin:
    output:
        source = "%s/custom/vancomycin.fasta" %database_dir
    log:
        stdout = "%s/Databases/vancomycin.log" %logdir
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
        source = "%s/custom/vancomycin_operon.fasta" %database_dir
    log:
        stdout = "%s/Databases/vancomycin_operon.log" %logdir
    message:
        "[fetch_vancomycin_operon]: Downloading custom vancomycin operon database"
    shell:
        """
        set -euo pipefail
        OUTDIR=$(dirname {output.source})
        mkdir -p $OUTDIR
                
        fasta_url="https://raw.githubusercontent.com/ssi-dk/ssi_analysis_utility_db/refs/heads/main/resistance/vancomycin_operon.fasta"

        cmd="curl -fSL $fasta_url -o {output.source}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval "$cmd" >> {log.stdout} 2>&1
        """