#------------------------------ Modules --------------------------------#

# Rule: PlasmidFinder
# Runs plasmidfinder analysis for a given sample using paired-end Illumina reads.
rule PlasmidFinder:
    input:
        # Input paired-end Illumina reads.
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        db_index = lambda wc: f"Logs/Databases/{species_configs[sample_to_organism[wc.sample]]['alignment_database']['plasmidfinder']['kma_index_flag']}.kma_index.done",
        db_dir = lambda wc: f"{getattr(rules, species_configs[sample_to_organism[wc.sample]]['alignment_database']['plasmidfinder']['db']).output.database}",
    output:
        # Output directory for plasmidfinder results.
        out_dir = directory(f"{OUT_FOLDER}" + "/{sample}/PlasmidFinder")
    conda:
        config["analysis_settings"]["plasmidfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/PlasmidFinder.log'
    message:
        "[PlasmidFinder]: Running PlasmidFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="plasmidfinder.py -i {input.R1} {input.R2} -o {output.out_dir} -p {input.db_dir} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# Rule: ResFinder
# Runs ResFinder to identify acquired resistance genes in the sample.
rule ResFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        res_db_index = lambda wc: f"Logs/Databases/{species_configs[sample_to_organism[wc.sample]]['alignment_database']['resfinder']['kma_index_flag']}.kma_index.done",
        res_db_dir = lambda wc: f"{getattr(rules, species_configs[sample_to_organism[wc.sample]]['alignment_database']['resfinder']['db']).output.database}",
        point_db_index = lambda wc: f"Logs/Databases/{species_configs[sample_to_organism[wc.sample]]['alignment_database']['pointfinder']['kma_index_flag']}.kma_index.done",
        point_db_dir = lambda wc: f"{getattr(rules, species_configs[sample_to_organism[wc.sample]]['alignment_database']['pointfinder']['db']).output.database}",
        disin_db_index = lambda wc: f"Logs/Databases/{species_configs[sample_to_organism[wc.sample]]['alignment_database']['disinfinder']['kma_index_flag']}.kma_index.done",
        disin_db_dir = lambda wc: f"{getattr(rules, species_configs[sample_to_organism[wc.sample]]['alignment_database']['disinfinder']['db']).output.database}",
    output:
        out_dir = directory(f"{OUT_FOLDER}" + "/{sample}/ResFinder")
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/ResFinder.log'
    message:
        "[ResFinder]: Running ResFinder, PointFinder, and DisinFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="python -m resfinder -ifq {input.R1} {input.R2} -o {output.out_dir} -db_res {input.res_db_dir} -db_disinf {input.disin_db_dir} -db_point {input.point_db_dir} -acq"
 
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# Rule: VirulenceFinder
# Runs VirulenceFinder to identify virulence genes in the sample.
rule VirulenceFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        vir_db_index = lambda wc: f"Logs/Databases/{species_configs[sample_to_organism[wc.sample]]['alignment_database']['virulencefinder']['kma_index_flag']}.kma_index.done",
        vir_db_dir = lambda wc: f"{getattr(rules, species_configs[sample_to_organism[wc.sample]]['alignment_database']['virulencefinder']['db']).output.database}",
    output:
        out_dir = directory(f"{OUT_FOLDER}" + "/{sample}/VirulenceFinder")
    conda:
        config["analysis_settings"]["virulencefinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/VirulenceFinder.log'
    message:
        "[VirulenceFinder]: Running VirulenceFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="virulencefinder.py -i {input.R1} {input.R2} -o {output.out_dir} -p {input.vir_db_dir}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


# Rule: SerotypeFinder
# Identifies serotypes from Illumina paired-end reads.
rule serotypefinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        ser_db_index = lambda wc: f"Logs/Databases/{species_configs[sample_to_organism[wc.sample]]['alignment_database']['serotypefinder']['kma_index_flag']}.kma_index.done",
        ser_db_dir = lambda wc: f"{getattr(rules, species_configs[sample_to_organism[wc.sample]]['alignment_database']['serotypefinder']['db']).output.database}",
        database = rules.setup_SerotypeFinder.output.database
    output:
        out_dir = directory(f"{OUT_FOLDER}" + "/{sample}/SerotypeFinder")
    conda:
        config["analysis_settings"]["serotypefinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/SerotypeFinder.log'
    message:
        "[SerotypeFinder]: Running SerotypeFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="serotypefinder -i {input.R1} {input.R2} -o {output.out_dir} -p {input.ser_db_dir} -x"

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
        result = "%s/{sample}/AMRFinder/AMR_{assembler}_{sample}.tsv" %OUT_FOLDER
    conda:
        config["analysis_settings"]["amrfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/AMRFinder/AMRFinder_{assembler}.log'
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