rule custom_kmeralignment:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_custom_kmeraligner_index.output.names
  params:
    prefix = "%s/{sample}/kmeraligner/{database}" %OUT_FOLDER
  output:
    results = "%s/{sample}/kmeraligner/{database}.res" %OUT_FOLDER,
    sam = temp("%s/{sample}/kmeraligner/{database}.sam" %OUT_FOLDER),
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

    db_path=$(dirname {input.database})/$(basename {input.database} .name)

    cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} -t_db $db_path -nc -nf -ref_fsa -sam 4 > {output.sam}"
    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


rule custom_kmeralignment_samtools_filtration:
  input:
    sam = rules.custom_kmeralignment.output.sam
  params:
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["samtools"]["view"]["options"]
  output:
    bam = temp("%s/{sample}/samtools/{database}.bam" %OUT_FOLDER)
  conda:
    config["analysis_settings"]["htslib"]["yaml"]
  log:
    stdout = "Logs/{sample}/custom_kmeralignment_samtools_filtration_{database}.log"
  message:
    "[custom_kmeralignment_samtools_filtration]: Filtering kmeralignment output for {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="samtools view {params.options} {input.sam} -o {output.bam}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


rule samtools_sort:
  input:
    bam = "%s/{sample}/samtools/{database}.bam" %OUT_FOLDER
  params:
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["samtools"]["sort"]["options"]
  output:
    bam_sort = temp("%s/{sample}/samtools/{database}_sorted.bam" %OUT_FOLDER)
  conda:
    config["analysis_settings"]["htslib"]["yaml"]
  log:
    stdout = "Logs/{sample}/samtools_sort_{database}.log"
  message:
    "[samtools_sort]: Sorting filtered bam for {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="samtools sort -o {output.bam_sort} {input.bam}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


rule samtools_index:
  input:
    bam_sort = rules.samtools_sort.output.bam_sort
  output:
    bam_index = temp("%s/{sample}/samtools/{database}_sorted.bam.bai" %OUT_FOLDER)
  conda:
    config["analysis_settings"]["htslib"]["yaml"]
  log:
    stdout = "Logs/{sample}/samtools_index_{database}.log"
  message:
    "[samtools_index]: Indexing sorted bam for {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="samtools index {input.bam_sort}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


rule bcftools_pileup:
  input:
    bam_sort = rules.samtools_sort.output.bam_sort,
    bam_index = rules.samtools_index.output.bam_index,
    reference = "%s/samtools/{database}.fasta" %database_path
  output:
    pileup = temp("%s/{sample}/bcftools/{database}.bcf" %OUT_FOLDER)
  conda:
    config["analysis_settings"]["htslib"]["yaml"]
  log:
    stdout = "Logs/{sample}/bcftools_pileup_{database}.log"
  message:
    "[bcftools_pileup]: Generating mpileup for {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="bcftools mpileup -Ob -f {input.reference} {input.bam_sort} -o {output.pileup}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """

rule bcftools_index:
  input:
    pileup = "%s/{sample}/bcftools/{database}.bcf" %OUT_FOLDER
  output:
    index = temp("%s/{sample}/bcftools/{database}.bcf.csi" %OUT_FOLDER)
  conda:
    config["analysis_settings"]["htslib"]["yaml"]
  log:
    stdout = "Logs/{sample}/bcftools_index_{database}.log"
  message:
    "[bcftools_index]: Indexing mpileup of {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="bcftools index -f {input.pileup}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


rule bcftools_filter_indels:
  input:
    pileup = rules.bcftools_pileup.output.pileup,
    index = rules.bcftools_index.output.index
  params:
    region = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["bcftools"]["view"]["region"],
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["bcftools"]["view"]["options"]
  output:
    indels = temp("%s/{sample}/bcftools/{database}_indels.bcf" %OUT_FOLDER)
  conda:
    config["analysis_settings"]["htslib"]["yaml"]
  log:
    stdout = "Logs/{sample}/bcftools_filter_indels_{database}.log"
  message:
    "[bcftools_filter_indels]: Filtering indels of {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="bcftools view -r {params.region} {params.options} -Ob -o {output.indels} {input.pileup}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


rule bcftools_variant_call:
  input:
    pileup = rules.bcftools_pileup.output.pileup,
    index = rules.bcftools_index.output.index
  output: 
    variants = temp("%s/{sample}/bcftools/{database}_variants.bcf" %OUT_FOLDER)
  conda:
    config["analysis_settings"]["htslib"]["yaml"]
  log:
    stdout = "Logs/{sample}/bcftools_variant_call_{database}.log"
  message:
    "[bcftools_variant_call]: Calling variant of {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="bcftools call -mv -Ob --ploidy 1 {input.pileup} -o {output.variants}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """