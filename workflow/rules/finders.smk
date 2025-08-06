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


# Rule: SerotypeFinder
# Identifies serotypes from Illumina paired-end reads.
rule serotypefinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_SerotypeFinder.output.database
    output:
        directory("%s/{sample}/SerotypeFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["serotypefinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/SerotypeFinder.log'
    message:
        "[SerotypeFinder]: Running SerotypeFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}

        cmd="serotypefinder -i {input.R1} {input.R2} -o {output} -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
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
