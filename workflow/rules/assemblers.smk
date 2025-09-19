#created a slight variation naming them according to the rules as well, such that it accomodates any potential seperate assembly rules
rule shovill:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
  output:
    assembly = "%s/shovill/{sample}/{assembler}.fasta" % OUT_FOLDER
  conda:
    "../envs/shovill.yaml"
  log:
    stdout = "Logs/shovill/{sample}_{assembler}.log"
  threads:
    lambda wc: max(1, int(workflow.cores * 0.66667))
  message:
    "[Shovill]: Assembling {wildcards.sample} using {wildcards.assembler} with {threads} cores. This may take some time!\nInspect {log.stdout} for more details!"
  shell:
    """
    outdir=$(dirname {output.assembly})
    assembler="$(echo "{wildcards.assembler}" | cut -d'_' -f2)"
    echo "shovill assembler: \n $assembler \n"

    cmd="shovill --R1 {input.R1} --R2 {input.R2} --outdir $outdir/$assembler --force --cpus {threads} --trim --assembler $assembler"
    echo "Executing command:\n$cmd\n" 

    echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1

    cmd_asm="rm $outdir/$assembler/$assembler.fasta && mv $outdir/$assembler/contigs.fa {output.assembly}"
    echo "Executing command:\n$cmd_asm\n"
    echo "Executing command:\n$cmd_asm\n" >> {log.stdout} 2>&1
    eval $cmd_asm >> {log.stdout} 2>&1
    """

#consider adding this:
#   params: extra = lambda wildcards: ("--noreadcorr --nocorr" if "skesa" in wildcards.assembler.lower() else "")
# also eventually change the output folder structure to:
#   assembly = "%s/{sample}/shovill/{assembler}.fasta" % OUT_FOLDER


# Rule: spades_assembly
rule spades:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    output:
        assembly = "%s/spades/{sample}/{assembler}.fasta" % OUT_FOLDER
    conda:
        "../envs/spades.yaml"
    log:
        stdout = "Logs/spades/{sample}_{assembler}.log"
    message:
      "[spades_assembly]: Perform assembly using spades on {wildcards.sample}, this will take some time!"
    threads:
      lambda wc: max(1, int(workflow.cores * 0.66667))
    shell:
      """
      outdir=$(dirname {output.assembly})
      assembler={wildcards.assembler}
      echo "assembler: \n $assembler \n"

      cmd="spades.py -1 {input.R1} -2 {input.R2} --threads {threads} --isolate --only-assembler -o $outdir"
      echo "Executing command:\n$cmd\n" 

      echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
      eval $cmd >> {log.stdout} 2>&1

      cmd_asm="mv $outdir/contigs.fasta {output.assembly}"
      echo "Executing command:\n$cmd_asm\n" >> {log.stdout} 2>&1
      eval $cmd_asm >> {log.stdout} 2>&1
      """

# Rule: Skesa_assembly
#  DeBruijn graph-based de-novo assembler for microbial genomes
rule skesa:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    output:
        assembly = "%s/skesa/{sample}/{assembler}.fasta" % OUT_FOLDER
    conda:
        "../envs/skesa.yaml"
    log:
        stdout = "Logs/skesa/{sample}_{assembler}.log"
    message:
        "[skesa_assembly]: Perform assembly using skesa on {wildcards.sample}, this will take some time!"
    threads:
      lambda wc: max(1, int(workflow.cores * 0.66667))
    shell:
      """
      cmd="skesa --reads {input.R1},{input.R2} --contigs_out {output.assembly} --cores {threads}"
      echo "Executing command:\n$cmd\n" 

      echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
      eval $cmd >> {log.stdout} 2>&1
      """