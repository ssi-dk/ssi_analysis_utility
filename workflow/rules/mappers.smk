rule custom_kmeralignment:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_custom_kmeraligner_index.output.names
  params:
    prefix_out = "%s/{sample}/kmeraligner/{database}" %output_folder,
    prefix_db = rules.setup_custom_kmeraligner_index.params.prefix    
  output:
    results = "%s/{sample}/kmeraligner/{database}.res" %output_folder,
    sam = temp("%s/{sample}/samtools/{database}.sam" %output_folder),
    done = touch(f"{output_folder}" + "/{sample}/kmeraligner/{database}.done")
  conda:
    "../envs/kmeraligner.yaml"
  log:
    stdout = "Logs/{sample}/custom_kmeralignment_{database}.log"
  message:
    "[kmeraligner]: Running KMA for {wildcards.database} on {wildcards.sample}"
  shell:
    """
    mkdir -p $(dirname {output.results})
    mkdir -p $(dirname {output.sam})

    cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix_out} -t_db {params.prefix_db} -na -nc -nf -sam 4 > {output.sam}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """

rule custom_kmerconsensus:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_custom_kmeraligner_index.output.names
  params:
    prefix_out = "%s/{sample}/kmerconsensus/{database}" %output_folder,
    prefix_db = rules.setup_custom_kmeraligner_index.params.prefix
  output:
    results = temp("%s/{sample}/kmerconsensus/{database}.res" %output_folder),
    seq = temp("%s/{sample}/kmerconsensus/{database}.fsa" %output_folder),
    aln = temp("%s/{sample}/kmerconsensus/{database}.aln" %output_folder)
  conda:
    "../envs/kmeraligner.yaml"
  log:
    stdout = "Logs/{sample}/custom_kmerconsensus_{database}.log"
  message:
    "[kmerconsensus]: Running KMA for {wildcards.database} on {wildcards.sample}"
  shell:
    """
    mkdir -p $(dirname {output.seq})

    cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix_out} -t_db {params.prefix_db} -nf -ref_fsa"
    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """

rule custom_bowtie2alignment:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_custom_bowtie2_index.output.bt2_1 # just locate one of the bt2 files to activate the db_setup
  params:
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["bowtie2"]["options"],
  output:
    sam = temp("%s/{sample}/bowtie2/{database}.sam" %output_folder)
  threads:
    max(1, workflow.cores * 0.3333333)
  conda:
    "../envs/bowtie2.yaml"
  log:
    stdout = "Logs/{sample}/custom_bowtie2_{database}.log"
  message:
    "[bowtie2aligner]: Running Bowtie2 for {wildcards.database} on {wildcards.sample} using {threads} thread(s)"
  shell:
    """
    mkdir -p $(dirname {output.sam})

    db_prefix="{input.database}"
    db_prefix="${{db_prefix%.1.bt2}}"
    
    cmd="bowtie2 -1 {input.R1} -2 {input.R2} -q -S {output.sam} {params.options} -x $db_prefix --threads {threads}"
    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """