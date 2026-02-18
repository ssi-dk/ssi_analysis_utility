rule spades:
    input:
        R1 = lambda wc: SAMPLESHEET.loc[wc.sample, "read1"],
        R2 = lambda wc: SAMPLESHEET.loc[wc.sample, "read2"]
    output:
        assembly = "%s/{sample}/spades/{sample}.fasta" %OUTDIR
    conda:
        "../envs/shovill.yaml"
    log:
        stdout = "%s/Assemblies/{sample}_spades.log" %LOGDIR
    threads:
        max(1, workflow.cores * 2 / 3)
    priority: 2
    message:
        "[SPAdes]: Assemblying {wildcards.sample} using SPAdes with {threads} thread(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        outdir=$(dirname {output.assembly})
        cmd="spades.py -1 {input.R1} -2 {input.R2} --threads {threads} --isolate -o $outdir"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="mv $outdir/contigs.fasta {output.assembly}"
        echo "### SPAdes Done!###\nExecuting command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule skesa:
    input:
        R1 = lambda wc: SAMPLESHEET.loc[wc.sample, "read1"],
        R2 = lambda wc: SAMPLESHEET.loc[wc.sample, "read2"]
    output:
        assembly = "%s/{sample}/skesa/{sample}.fasta" %OUTDIR
    conda:
        "../envs/shovill.yaml"
    log:
        stdout = "%s/Assemblies/{sample}_Skesa.log" %LOGDIR
    threads:
        max(1, workflow.cores * 2 / 3)
    priority: 2
    message:
        "[Skesa]: Assemblying {wildcards.sample} using Skesa with {threads} core(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        cmd="skesa --reads {input.R1},{input.R2} --contigs_out {output.assembly} --cores {threads}"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule shovill:
    input:
        R1 = lambda wc: SAMPLESHEET.loc[wc.sample, "read1"],
        R2 = lambda wc: SAMPLESHEET.loc[wc.sample, "read2"]
    output:
        assembly = "%s/{sample}/shovill/{sample}.fasta" %OUTDIR
    conda:
        "../envs/shovill.yaml"
    log:
        stdout = "%s/Assemblies/{sample}_Shovill.log" %LOGDIR
    threads:
        max(1, workflow.cores * 2  / 3)
    priority: 2
    message:
        "[Shovill]: Assemblying {wildcards.sample} using Shovill with {threads} CPU(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        mkdir -p $(dirname {output.assembly})
        outdir=$(dirname {output.assembly})

        cmd="shovill --R1 {input.R1} --R2 {input.R2} --outdir $outdir/ --force --cpus {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="mv $outdir/contigs.fa {output.assembly}"
        echo "### Shovil Done!###\nExecuting command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule assembly:
    input:
        input_assembly = "%s/{sample}/{assembler}/{sample}.fasta" %OUTDIR
    output:
        output_assembly = "%s/{sample}/Assemblies/{sample}_{assembler}.fasta" %OUTDIR
    log:
        stdout = "%s/Assemblies/{sample}_{assembler}_assembly.log" %LOGDIR
    shell:
        """
        mkdir -p $(dirname {output.output_assembly})

        cmd="cp {input.input_assembly} {output.output_assembly}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 
        """
