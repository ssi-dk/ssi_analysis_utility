rule spades:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = f"{outdir}/{{sample}}/raw/spades/spades_{{sample}}.fasta"
    conda:
        ENVS_DIR / "spades.yaml"
    log:
        stdout = f"{logdir}/spades_{{sample}}.log"
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 2
    message:
        "[spades]: Assemblying {wildcards.sample} using SPAdes with {threads} thread(s). This may take some time!\nInspect {log.stdout} for more details!"
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
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = f"{outdir}/{{sample}}/raw/skesa/skesa_{{sample}}.fasta"
    conda:
        ENVS_DIR / "skesa.yaml"
    log:
        stdout = f"{logdir}/skesa_{{sample}}.log"
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 2
    message:
        "[skesa]: Assemblying {wildcards.sample} using Skesa with {threads} core(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        cmd="skesa --reads {input.R1},{input.R2} --contigs_out {output.assembly} --cores {threads}"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule shovill:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = f"{outdir}/{{sample}}/raw/shovill/shovill_{{sample}}.fasta"
    conda:
        ENVS_DIR / "shovill.yaml"
    log:
        stdout = f"{logdir}/shovill_{{sample}}.log"
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 2
    message:
        "[shovill]: Assemblying {wildcards.sample} using Shovill with {threads} CPU(s). This may take some time!\nInspect {log.stdout} for more details!"
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
        assembly = f"{outdir}/{{sample}}/raw/{{assembler}}/{{assembler}}_{{sample}}.fasta"
    output:
        assembly = f"{outdir}/{{sample}}/Assemblies/{{assembler}}_{{sample}}.fasta"
    log:
        stdout = f"{logdir}/sync_{{assembler}}_{{sample}}.log"
    message:
        "[assembly]: Syncronizing {wildcards.assembler} for {wildcards.sample} from raw location to assembly folder"
 
    shell:
        """
        mkdir -p $(dirname {output.assembly})

        cmd="rsync -av {input.assembly} {output.assembly}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 
        """
