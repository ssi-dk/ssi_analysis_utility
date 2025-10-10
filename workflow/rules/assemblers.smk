rule spades:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
  output:
    assembly = "%s/{sample}/spades/{sample}.fasta" %output_folder,
  conda:
    "../envs/shovill.yaml"
  log:
    stdout = "Logs/Assemblies/{sample}_spades.log"
  threads:
    max(1, workflow.cores * 0.66667)
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
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
  output:
    assembly = "%s/{sample}/skesa/{sample}.fasta" %output_folder,
  conda:
    "../envs/shovill.yaml"
  log:
    stdout = "Logs/Assemblies/{sample}_Skesa.log"
  threads:
    max(1, workflow.cores * 0.66667)
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
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
  output:
    assembly = "%s/{sample}/shovill/{sample}.fasta" %output_folder,
  conda:
    "../envs/shovill.yaml"
  log:
    stdout = "Logs/Assemblies/{sample}_Shovill.log"
  threads:
    max(1, workflow.cores * 0.66667)
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
    input_assembly = "%s/{sample}/{assembler}/{sample}.fasta" %output_folder
  output:
    output_assembly = "%s/Assemblies/{sample}_{assembler}.fasta" %output_folder,
    done = temp("%s/Assemblies/{sample}_{assembler}.done" % output_folder)
  log:
    stdout = "Logs/Assemblies/{sample}_{assembler}_assembly.log"
  shell:
    """
    mkdir -p $(dirname {output.output_assembly})

    cmd="cp {input.input_assembly} {output.output_assembly}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1 

    echo "Assembly successfully created for {wildcards.sample}" > {output.done}
    """
