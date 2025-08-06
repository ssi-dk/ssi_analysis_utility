# Rule MLST:
# Runs Multi Locus Sequence Type to determine the ST profile of isolate
rule MLST:
    input:
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample],
        datefile = rules.update_MLST.output.datefile
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
