#------------------------------ Modules --------------------------------#

# Rule: PlasmidFinder
# Runs plasmidfinder analysis for a given sample using paired-end Illumina reads.
rule PlasmidFinder:
    input:
        # Input paired-end Illumina reads.
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_PlasmidFinder.output.database
    output:
        # Output directory for plasmidfinder results.
        directory("%s/{sample}/PlasmidFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["plasmidfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/PlasmidFinder.log'
    message:
        "[PlasmidFinder]: Running PlasmidFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}

        cmd="plasmidfinder.py -i {input.R1} {input.R2} -o {output} -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# Rule: ResFinder
# Runs ResFinder to identify acquired resistance genes in the sample.
rule ResFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        res_database = rules.setup_ResFinder.output.database,
        point_database = rules.setup_PointFinder.output.database,
        disin_database = rules.setup_DisinFinder.output.database
    output:
        directory("%s/{sample}/ResFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/ResFinder.log'
    message:
        "[ResFinder]: Running ResFinder, PointFinder, and DisinFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}

        cmd="python -m resfinder -ifq {input.R1} {input.R2} -o {output} -db_res {input.res_database} -db_disinf {input.disin_database} -db_point {input.point_database} -acq"
 
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# Rule: VirulenceFinder
# Runs VirulenceFinder to identify virulence genes in the sample.
rule VirulenceFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_VirulenceFinder.output.database
    output:
        directory("%s/{sample}/VirulenceFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["virulencefinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/VirulenceFinder.log'
    message:
        "[VirulenceFinder]: Running VirulenceFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}

        cmd="virulencefinder.py -i {input.R1} {input.R2} -o {output} -p {input.database}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# Rule: LREFinder
# Runs LRE-Finder to analyze long-read sequencing data.
rule LREFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the LRE database, application, and additional options.
        db_path = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["LRE-finder"]["database"],
        app_path = config["analysis_settings"]["LRE-finder"]["script"],
        min_con_ID = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["LRE-finder"]["min_consensus_ID"],
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["LRE-finder"]["additional_option"],
        prefix = "OUT_FOLDER/{sample}/lre-finder/{sample}"
    conda:
        config["analysis_settings"]["LRE-finder"]["yaml"]
    output:
        directory("OUT_FOLDER/{sample}/lre-finder/")
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
        
        # Run LRE-Finder with the specified parameters and inputs.
        python {params.app_path}/LRE-Finder.py -ipe {input.R1} {input.R2} -o {params.prefix} -t_db {params.db_path} -ID {params.min_con_ID} {params.add_opt} 
        """

# Rule: serotypefinder
# Identifies serotypes from Illumina paired-end reads.
rule serotypefinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the serotype database.
        db_path = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["serotypefinder"]["database"],
    output:
        directory("OUT_FOLDER/{sample}/serotypefinder/")
    conda:
        config["analysis_settings"]["serotypefinder"]["yaml"]
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
# rule kmerfinder:
#     input:
#         R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
#         R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
#     params:
#         # Path to the kmerfinder database, KMA aligner, and taxa file.
#         db_path = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kmerfinder"]["database"],
#         taxa = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kmerfinder"]["taxa"]
#     output:
#         directory("OUT_FOLDER/{sample}/kmerfinder/")
#     conda:
#         config["analysis_settings"]["kmerfinder"]["yaml"]
#     shell:
#         """
#         # Check if the output directory exists, and skip if it does.
#         if [ -d {output} ]; 
#             then
#                 echo "Directory {output} exists, skipping."
#                 exit 1
#             else
#                 mkdir {output}
#         fi 
#         # Run kmerfinder.py with the specified parameters and inputs.
#         kmerfinder.py  -i {input.R1} {input.R2} -o {output} -db {params.db_path} -tax {params.taxa}
#         """

# Rule: cgMLSTFinder
# Runs cgMLSTFinder for comparative genomic typing.
rule cgMLSTFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the cgMLSTFinder app, database, and KMA aligner.
        app_path = config["analysis_settings"]["cgMLSTFinder"]["script"],
        db_path = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["cgMLSTFinder"]["database"],
        scheme = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["cgMLSTFinder"]["scheme"],
    output:
        directory("OUT_FOLDER/{sample}/cgmlstfinder/")
    conda:
        config["analysis_settings"]["cgMLSTFinder"]["yaml"]   
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
        python {params.app_path}/cgMLST.py  -i {input.R1},{input.R2} -o {output} -db {params.db_path}  -s {params.scheme}
        """

# Rule: AMRFinder
# Runs AMRFinder to identify acquired resistance genes in the sample.
rule AMRFinder:
    input:
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample],
        database = rules.setup_AMRFinder.output.database
    params:
        # Point mutation
        organism = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["amrfinder"]["organism"]
    output:
        # "{out}/{sample}/amrfinder/{sample}.tsv"
        directory("%s/{sample}/AMRFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["amrfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/AMRFinder.log'
    message:
        "[AMRFinder]: Running AMRFinderFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}

        cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.organism} --output {output}/{wildcards.sample}.tsv"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """
