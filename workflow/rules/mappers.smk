rule custom_kmeralignment:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_custom_kmeraligner_index.output.names
  params:
    prefix = "%s/{sample}/kmeraligner/{database}" %OUT_FOLDER
  output:
    results = "%s/{sample}/kmeraligner/{database}.res" %OUT_FOLDER,
    sam = temp("%s/{sample}/samtools/{database}.sam" %OUT_FOLDER),
    seq = temp("%s/{sample}/kmeraligner/{database}.fsa" %OUT_FOLDER)
  conda:
    config["analysis_settings"]["KMA"]["yaml"]
  log:
    stdout = "Logs/{sample}/custom_kmeralignment_{database}.log"
  message:
    "[kmeraligner]: Running KMA for {wildcards.database} on {wildcards.sample}"
  shell:
    """
    mkdir -p $(dirname {output.results})
    mkdir -p $(dirname {output.sam})
    

    db_path=$(dirname {input.database})/$(basename {input.database} .name)

    cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} -t_db $db_path -nc -nf -ref_fsa -sam 4 > {output.sam}"
    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


