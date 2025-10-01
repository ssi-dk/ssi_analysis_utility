rule shovill:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
  output:
    assembly = "%s/Assemblies/{sample}/{assembler}.fasta" %output_folder
  conda:
    "../envs/shovill.yaml"
  log:
    stdout = "Logs/Assemblies/{sample}_{assembler}.log"
  threads:
    workflow.cores * 0.66667
  message:
    "[Shovill]: Assemblying {wildcards.sample} using {wildcards.assembler} with {threads} cores. This may take some time!\nInspect {log.stdout} for more details!"
  shell:
    """
    outdir=$(dirname {output.assembly})
    assembler={wildcards.assembler}

    cmd="shovill --R1 {input.R1} --R2 {input.R2} --outdir $outdir/$assembler --force --cpus {threads} --assembler $assembler --trim && mv $outdir/$assembler/contigs.fa {output.assembly}"

    echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """
