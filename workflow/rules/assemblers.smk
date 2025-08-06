# Rule SPAdes
# Runs SPAdes to assemble reads into cotigs
rule spades:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]

    output:
        assembly = "%s/{sample}/spades/{sample}.fasta" %OUT_FOLDER

    conda:
        config["analysis_settings"]["spades"]["yaml"]

    log:
        stdout = "Logs/{sample}/spades.log"
    
    message:
        "[SPAdes]: Running SPAdes on {wildcards.sample}"

    threads:
        workflow.cores/4
    shell:
        """
        cmd="spades.py -1 {input.R1} -2 {input.R2} --isolate --threads {threads} -o {output.assembly}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 

        cmd="mv $(dirname {output.assembly})/contigs.fasta $(dirname {output.assembly})/{wildcards.sample}.fasta"
        echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1         
        """

