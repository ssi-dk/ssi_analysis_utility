rule setup_PlasmidFinder:
    conda:
        config["analysis_settings"]["plasmidfinder"]["yaml"]
    output: 
        database = directory(f'{database_path}/{config["analysis_settings"]["plasmidfinder"]["database"]}')
    log:
        stdout = f'Logs/plasmidfinder_db.log'
    message:
        "[PlasmidFinder_db]: Setting up PlasmidFinder_db"
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
        stdout = f'Logs/resfinder_db.log'
    message:
        "[ResFinder_db]: Setting up ResFinder_db"
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
        stdout = f'Logs/pointfinder_db.log'
    message:
        "[PointFinder_db]: Setting up PointFinder_db"
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
        stdout = f'Logs/disinfinder_db.log'
    message:
        "[DisinFinder_db]: Setting up DisinFinder_db"
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
        stdout = f'Logs/virulencefinder_db.log'
    message:
        "[VirulenceFinder_db]: Setting up VirulenceFinder_db"
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
