rule custom_kmeralignment:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_custom_kmeraligner_index.output.names
  params:
    prefix = "%s/{sample}/kmeraligner/{database}" %OUT_FOLDER
  output:
    results = "%s/{sample}/kmeraligner/{database}.res" %OUT_FOLDER,
    sam = temp("%s/{sample}/samtools/{database}.sam" %OUT_FOLDER)
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

    db_path=$(dirname {input.database})/$(basename {input.database} .name)

    cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} -t_db $db_path -na -nc -nf -sam 4 > {output.sam}"
    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """

rule custom_kmerconsensus:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_custom_kmeraligner_index.output.names
  params:
    prefix = "%s/{sample}/kmerconsensus/{database}" %OUT_FOLDER
  output:
    results = temp("%s/{sample}/kmerconsensus/{database}.res" %OUT_FOLDER),
    seq = "%s/{sample}/kmerconsensus/{database}.fsa" %OUT_FOLDER
  conda:
    "../envs/kmeraligner.yaml"
  log:
    stdout = "Logs/{sample}/custom_kmerconsensus_{database}.log"
  message:
    "[kmerconsensus]: Running KMA for {wildcards.database} on {wildcards.sample}"
  shell:
    """
    mkdir -p $(dirname {output.seq})   

    db_path=$(dirname {input.database})/$(basename {input.database} .name)

    cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} -t_db $db_path -nf -ref_fsa"
    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """

rule custom_bowtie2alignment:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_custom_bowtie2aligner_index.output.bt2_1 # just locate one of the bt2 files to activate the db_setup
  params:
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["bowtie2aligner"]["options"],
  output:
    sam = "%s/{sample}/bowtie2aligner/{database}.sam" %OUT_FOLDER,
  threads:
    min(workflow.cores, 16)
  conda:
    "../envs/bowtie2aligner.yaml"
  log:
    stdout = "Logs/{sample}/custom_bowtie2_{database}.log"
  message:
    "[bowtie2aligner]: Running Bowtie2 for {wildcards.database} on {wildcards.sample}"
  shell:
    """
    mkdir -p $(dirname {output.sam})

    db_prefix="{input.database}"
    db_prefix="${{db_prefix%.1.bt2}}"
    
    cmd="bowtie2 -1 {input.R1} -2 {input.R2} -q -S {output.sam} {params.options} -x $db_prefix --threads {threads}"
    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """