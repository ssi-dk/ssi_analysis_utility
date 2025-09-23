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
    seq = temp("%s/{sample}/kmerconsensus/{database}.fsa" %OUT_FOLDER),
    aln = temp("%s/{sample}/kmerconsensus/{database}.aln" %OUT_FOLDER)
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
    database = rules.setup_custom_bowtie2_index.output.bt2_1 # just locate one of the bt2 files to activate the db_setup
  params:
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["bowtie2"]["options"],
  output:
    sam = temp("%s/{sample}/bowtie2/{database}.sam" %OUT_FOLDER)
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