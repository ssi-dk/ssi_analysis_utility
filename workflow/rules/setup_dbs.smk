rule setup_PlasmidFinder:
    conda:
        "../envs/plasmidfinder.yaml"
    output: 
        database = directory("%s/plasmidfinder_db" %database_path),
        version_db = "%s/plasmidfinder_db/PlasmidFinder_version.txt" %database_path
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
    conda:
        "../envs/resfinder.yaml"
    output:
        database = directory("%s/resfinder_db" %database_path),
        version_db = "%s/resfinder_db/ResFinder_version.txt" %database_path
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
    conda:
        "../envs/resfinder.yaml"
    output:
        database = directory("%s/pointfinder_db" %database_path),
        version_db = "%s/pointfinder_db/PointFinder_version.txt" %database_path
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
    conda:
        "../envs/resfinder.yaml"
    output:
        database = directory("%s/disinfinder_db" %database_path),
        version_db = "%s/disinfinder_db/DisinFinder_version.txt" %database_path
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
    conda:
        "../envs/virulencefinder.yaml"
    output:
        database = directory("%s/virulencefinder_db" %database_path),
        version_db = "%s/virulencefinder_db/VirulenceFinder_version.txt" %database_path
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
    conda:
        "../envs/serotypefinder.yaml"
    output:
        database = directory("%s/serotypefinder_db" %database_path),
        version_db = "%s/serotypefinder_db/SerotypeFinder_version.txt" %database_path
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
        database = directory("%s/spatyper_db" %database_path)
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
    conda:
        "../envs/amrfinder.yaml"
    output:
        database = directory("%s/amrfinderplus/latest" %database_path),
        version_db = "%s/amrfinderplus/latest/AMRFinder_version.txt" %database_path
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
        version_cmd="amrfinder_update -v"
        date_cmd="date -I"
        
        echo -e "Executing command:\n$modified_cmd\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        modified_str="$(eval "$modified_cmd" 2>> {log.stdout})"
        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s%s%s\t%s\n' "AMRFinder_" "$version_str" "_" "$modified_str" "$date_str" > {output.version_db}
        """


rule setup_custom_kmeraligner_index:
  input:
    source = "%s/custom/{database}.fasta" %database_path
  params:
    prefix = "%s/kmeraligner/{database}" %database_path
  output:
    combined_size = "%s/kmeraligner/{database}.comp.b" %database_path,
    lengths = "%s/kmeraligner/{database}.length.b" %database_path,
    names = "%s/kmeraligner/{database}.name" %database_path,
    seqs = "%s/kmeraligner/{database}.seq.b" %database_path,
    version_db = "%s/kmeraligner/{database}_kmaindex_version.txt" %database_path
  conda:
    "../envs/kmeraligner.yaml"
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
    source = "%s/custom/{database}.fasta" %database_path
  params:
    prefix = "%s/bowtie2/{database}" %database_path
  output:
    bt2_1 = "%s/bowtie2/{database}.1.bt2" %database_path,
    bt2_2 = "%s/bowtie2/{database}.2.bt2" %database_path,
    bt2_3 = "%s/bowtie2/{database}.3.bt2" %database_path,
    bt2_4 = "%s/bowtie2/{database}.4.bt2" %database_path,
    bt2_1_rev = "%s/bowtie2/{database}.rev.1.bt2" %database_path,
    bt2_2_rev = "%s/bowtie2/{database}.rev.2.bt2" %database_path,
    version_db = "%s/bowtie2/{database}_bowtie2index_version.txt" %database_path
  conda:
    "../envs/bowtie2.yaml"
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
    source = "%s/custom/{database}.fasta" %database_path
  output:
    source = "%s/samtools/{database}.fasta" %database_path, 
    index = "%s/samtools/{database}.fasta.fai" %database_path
  conda:
    "../envs/htslib.yaml"
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

# rule update_databases:
#     input:
#         rules.setup_AMRFinder.output,
#         rules.setup_custom_bowtie2_index.output,
#         rules.setup_custom_kmeraligner_index.output,
#         rules.setup_custom_samtool_index.output,
#         rules.setup_DisinFinder.output,
#         rules.setup_PlasmidFinder.output,
#         rules.setup_PointFinder.output,
#         rules.setup_ResFinder.output,
#         rules.setup_SerotypeFinder.output,
#         rules.setup_VirulenceFinder.output
