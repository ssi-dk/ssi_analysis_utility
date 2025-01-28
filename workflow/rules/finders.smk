#------------------------------ Modules --------------------------------#

# Rule: PlasmidFinder
# Runs plasmidfinder analysis for a given sample using paired-end Illumina reads.
rule PlasmidFinder:
    input:
        # Input paired-end Illumina reads.
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the plasmid database and KMA aligner.
        db_path = config["analysis_settings"]["plasmidfinder"]["database"],
        kma_path = config["analysis_settings"]["kma_aligner"]["path"]
    output:
        # Output directory for plasmidfinder results.
        directory("{out}/{sample}/plasmidfinder/")
    shell:
        """
        # Check if the output directory exists, and skip if it does.
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1
            else
                mkdir {output}
        fi        
        # Run plasmidfinder.py with specified input, output, and parameters.
        plasmidfinder.py -i {input.R1} {input.R2} -o {output} -p {params.db_path} -mp {params.kma_path} -x
        """

# Rule: ResFinder
# Runs ResFinder to identify acquired resistance genes in the sample.
rule ResFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Paths to the resistance, point mutation, and disinfectant resistance databases.
        res_db_path = config["analysis_settings"]["resfinder"]["resfinder_db"],
        point_db_path = config["analysis_settings"]["resfinder"]["pointfinder_db"],
        disi_db_path = config["analysis_settings"]["resfinder"]["disinfinder_db"]
    output:
        directory("{out}/{sample}/resfinder/")
    shell:
        """
        # Check if the output directory exists, and skip if it does.
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1
            else
                mkdir {output}
        fi        
        # Run resfinder.py with specified input, output, and parameters.
        python -m resfinder -ifq {input.R1} {input.R2} -o {output} -db_res {params.res_db_path} -db_disinf {params.disi_db_path} -db_point {params.point_db_path} -acq
        """

# Rule: VirulenceFinder
# Runs VirulenceFinder to identify virulence genes in the sample.
rule VirulenceFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the virulence gene database and KMA aligner.
        db_path = config["analysis_settings"]["virulencefinder"]["database"],
        kma_path = config["analysis_settings"]["kma_aligner"]["path"]
    output:
        directory("{out}/{sample}/virulencefinder/")
    shell:
        """
        # Check if the output directory exists, and skip if it does.
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1
            else
                mkdir {output}
        fi 

        # Run virulencefinder.py with specified input, output, and parameters.
        virulencefinder.py -i {input.R1} {input.R2} -o {output} -p {params.db_path} -mp {params.kma_path}
        """

# Rule: LREFinder
# Runs LRE-Finder to analyze long-read sequencing data.
rule LREFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the LRE database, application, and additional options.
        db_path = config["analysis_settings"]["LRE-finder"]["database"],
        app_path = config["analysis_settings"]["LRE-finder"]["path"],
        min_con_ID = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["lre-finder"]["min_consensus_ID"],
        add_opt = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["lre-finder"]["additional_option"]
    output:
        directory("{out}/{sample}/lre-finder/")
    shell:
        """
        # Check if the output directory exists, and skip if it does.
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1
            else
                mkdir {output}
                cd {output}
        fi 
        # Run LRE-Finder with the specified parameters and inputs.
        python {params.app_path}/LRE-Finder.py -ipe {input.R1} {input.R2} -o lre-finder -t_db {params.db_path} -ID {params.min_con_ID} {params.add_opt} 
        """

# Rule: serotypefinder
# Identifies serotypes from Illumina paired-end reads.
rule serotypefinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the serotype database.
        db_path = config["analysis_settings"]["serotypefinder"]["database"],
    output:
        directory("{out}/{sample}/serotypefinder/")
    shell:
        """
        # Check if the output directory exists, and skip if it does.
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1
            else
                mkdir {output}
        fi 
        # Run serotypefinder with the specified inputs and database.
        serotypefinder -i {input.R1} {input.R2} -o {output} -p {params.db_path} 
        """

# Rule: kmerfinder
# Identifies microbial species or strain using k-mer-based alignment.
rule kmerfinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the kmerfinder database, KMA aligner, and taxa file.
        db_path = config["analysis_settings"]["kmerfinder"]["database"],
        kma_path = config["analysis_settings"]["kma_aligner"]["path"],
        taxa = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["kmerfinder"]["taxa"]
    output:
        directory("{out}/{sample}/kmerfinder/")
    shell:
        """
        # Check if the output directory exists, and skip if it does.
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1
            else
                mkdir {output}
        fi 
        # Run kmerfinder.py with the specified parameters and inputs.
        kmerfinder.py  -i {input.R1} {input.R2} -o {output} -db {params.db_path} -tax {params.taxa}
        """

# Rule: cgMLSTFinder
# Runs cgMLSTFinder for comparative genomic typing.
rule cgMLSTFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the cgMLSTFinder app, database, and KMA aligner.
        app_path = config["analysis_settings"]["cgMLSTFinder"]["path"],
        db_path = config["analysis_settings"]["cgMLSTFinder"]["database"],
        scheme = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["cgMLSTFinder"]["scheme"],
        kma_path = config["analysis_settings"]["kma_aligner"]["path"]
    output:
        directory("{out}/{sample}/cgmlstfinder/")
    shell:
        """
        # Check if the output directory exists, and skip if it does.
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1
            else
                mkdir {output}
        fi 
        # Run cgMLSTFinder with the specified inputs, outputs, and scheme.
        python {params.app_path}/cgMLST.py  -i {input.R1},{input.R2} -o {output} -db {params.db_path} -k {params.kma_path} -s {params.scheme}
        """
