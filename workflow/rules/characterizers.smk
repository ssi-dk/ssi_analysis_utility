# Rule MLST:
# Runs Multi Locus Sequence Type to determine the ST profile of isolate
rule MLST:
    input:
        assembly = rules.shovill.output.assembly,
        datefile = rules.update_MLST.output.datefile
    output:
        # "%s/{sample}/MLST/{sample}.tsv" %OUT_FOLDER
        mlst_file = "%s/{sample}/MLST/{assembler}_mlst.tsv" %OUT_FOLDER
    conda:
        "../envs/mlst.yaml"
    log:
    	stdout = "Logs/{sample}/{assembler}_mlst.log"
    message:
    	"[MLST]: Running MLST on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.mlst_file})

        cmd="mlst {input.assembly} --label $(basename {input.assembly} .fasta) > {output.mlst_file}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    	"""

# Rule Kleborate:
# Runs Kleborate characterising virulence and resistance in pathogen assemblies
rule kleborate:
    input:
        assembly = rules.shovill.output.assembly
    output:
       kleborate_outdir = directory("%s/{sample}/Kleborate/{assembler}" %OUT_FOLDER)
    params:
        options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kleborate"]["options"]
    conda:
        config["analysis_settings"]["kleborate"]["yaml"]
    log:
    	stdout = "Logs/{sample}/Kleborate_{assembler}.log"
    message:
    	"Kleborate: Running Kleborate on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        echo {config}
        mkdir -p $(dirname {output.kleborate_outdir})

        cmd="kleborate -a {input.assembly} --outdir {output.kleborate_outdir} {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    	"""

# Rule meningotype
# Runs meningotype on SPAdes assembled contigs to perform serotyping of N. Meningmeningitidis contigs
rule meningotype:
    input:
        assembly = rules.shovill.output.assembly
    output:
        meningotype = "%s/{sample}/meningotype/{assembler}_meningotype.tsv" %OUT_FOLDER
    conda:
        "../envs/meningotype.yaml"
    log:
        stdout = "Logs/{sample}/{assembler}_meningotype.log"
    message:
    	"[Meningotype]: Running Meningotype on {wildcards.assembler} assembly for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.meningotype})

        cmd="meningotype --all {input.assembly} > {output.meningotype} 2> {log.stdout}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    	"""

