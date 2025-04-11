rule setup_PlasmidFinder:
    conda:
        config["analysis_settings"]["plasmidfinder"]["yaml"]
    output: 
        database = directory(f'{database_path}/{config["analysis_settings"]["plasmidfinder"]["database"]}')
    log:
        stdout = f'Logs/setup_PlasmidFinder.log'
    message:
        "[setup_PlasmidFinder]: Setting up PlasmidFinder database"
    shell:
        """
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
        """


rule setup_ResFinder:
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["resfinder"]["database"]}')
    log:
        stdout = f'Logs/setup_ResFinder.log'
    message:
        "[setup_ResFinder]: Setting up ResFinder database"
    shell:
        """
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
        """


rule setup_PointFinder:
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["pointfinder"]["database"]}')
    log:
        stdout = f'Logs/setup_PointFinder.log'
    message:
        "[setup_PointFinder]: Setting up PointFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1

            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[pointfinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the pointfinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        """


rule setup_DisinFinder:
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["disinfinder"]["database"]}')
    log:
        stdout = f'Logs/setup_DisinFinder.log'
    message:
        "[setup_DisinFinder]: Setting up DisinFinder database"
    shell:
        """
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
        """


rule setup_VirulenceFinder:
    conda:
        config["analysis_settings"]["virulencefinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["virulencefinder"]["database"]}')
    log:
        stdout = f'Logs/setup_VirulenceFinder.log'
    message:
        "[setup_VirulenceFinder]: Setting up VirulenceFinder database"
    shell:
        """
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
        """


rule setup_AMRFinder:
    conda:
        config["analysis_settings"]["amrfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["amrfinder"]["database"]}')
    log:
        stdout = f'Logs/setup_AMRFinder.log'
    message:
        "[setup_AMRFinder]: Setting up AMRFinderPlus database"
    shell:
        """
        amrfinder_update --database $(dirname {output.database})
        """
