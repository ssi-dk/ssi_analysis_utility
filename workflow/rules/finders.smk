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
        out_dir = directory("%s/{sample}/PlasmidFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["plasmidfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/PlasmidFinder.log'
    message:
        "[PlasmidFinder]: Running PlasmidFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="plasmidfinder.py -i {input.R1} {input.R2} -o {output.out_dir} -p {input.database} -x"

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
        #point_database = rules.setup_PointFinder.output.database, #Pointfinder requires `species` definition
        disin_database = rules.setup_DisinFinder.output.database
    output:
        out_dir = directory("%s/{sample}/ResFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/ResFinder.log'
    message:
        "[ResFinder]: Running ResFinder, PointFinder, and DisinFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o {output} --acquired -db_res {input.res_database} --disinfectant -db_disinf {input.disin_database}"
 
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
        out_dir = directory("%s/{sample}/VirulenceFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["virulencefinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/VirulenceFinder.log'
    message:
        "[VirulenceFinder]: Running VirulenceFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="virulencefinder.py -i {input.R1} {input.R2} -o {output.out_dir} -p {input.database}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


# Rule: SerotypeFinder
# Identifies serotypes from Illumina paired-end reads.
rule serotypefinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_SerotypeFinder.output.database
    output:
        out_dir = directory("%s/{sample}/SerotypeFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["serotypefinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/SerotypeFinder.log'
    message:
        "[SerotypeFinder]: Running SerotypeFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="serotypefinder -i {input.R1} {input.R2} -o {output.out_dir} -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# Rule: AMRFinder
# Runs AMRFinder to identify acquired resistance genes in the sample.
rule AMRFinder:
    input:
        assembly = lambda wildcards: os.path.join(
            OUT_FOLDER,
            wildcards.sample,
            wildcards.assembler,
            {
                "spades": "contigs.fasta",
                "skesa": f"{wildcards.sample}.contigs.fasta"
            }[wildcards.assembler]
        ),
        database = rules.setup_AMRFinder.output.database
    params:
        # Point mutation
        organism = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["amrfinder"]["organism"]
    output:
        result = "%s/{sample}/AMRFinder/{assembler}.tsv" %OUT_FOLDER
    conda:
        config["analysis_settings"]["amrfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/AMRFinder_{assembler}.log'
    message:
        "[AMRFinder]: Running AMRFinderFinder for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        mkdir -p $(dirname {output.result})

        cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.organism} --output {output.result}"

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
        directory("%s/{sample}/lre-finder/" %OUT_FOLDER)
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
