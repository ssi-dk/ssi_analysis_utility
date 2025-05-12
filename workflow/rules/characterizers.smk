# Rule MLST:
# Runs Multi Locus Sequence Type to determine the ST profile of isolate
rule MLST:
    input:
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample],
    output:
        # mlst_file = "%s/{sample}/MLST/{sample}.tsv" %OUT_FOLDER
        directory("%s/{sample}/MLST" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["mlst"]["yaml"]
    log:
    	stdout = "Logs/{sample}/MLST.log"
    message:
    	"[MLST]: Running MLST on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}

        cmd="mlst {input.assembly} --label $(basename {input.assembly} .fasta) > {output}/{wildcards.sample}.tsv"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    	"""

# Rule Legionella_SBT:
# Runs Multi Locus Sequence Type to determine the ST profile of isolate
rule Legionella_SBT:
    input:
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample],
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]

    output:
        directory("%s/{sample}/Legionella_SBT" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["Legionella_SBT"]["yaml"]
    log:
    	stdout = "Logs/{sample}/Legionella_SBT.log"
    message:
    	"[Legionella SBT]: Running Legionella SBT on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}

        cmd="legionella_sbt --assembly_file {input.assembly} --r1_file {input.R1} --r2_file {input.R2} --output_folder {output}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    	"""
       

# Rule meningotype:
# Runs Multi Locus Sequence Type to determine the ST profile of isolate
rule meningotype:
    input:
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample],
    output:
        directory("%s/{sample}/meningotype" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["meningotype"]["yaml"]
    log:
    	stdout = "Logs/{sample}/meningotype.log"
    message:
    	"[Meningotype]: Running Meningotype on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}

        cmd="meningotype --all {input.assembly} > {output}/{wildcards.sample}.tsv"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    	"""

        