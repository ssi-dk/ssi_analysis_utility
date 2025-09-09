rule custom_kmeralignment:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_custom_kmeraligner_index.output.names
  params:
    prefix = "%s/{sample}/kmeraligner/{database}" %OUT_FOLDER
  output:
    results = "%s/{sample}/kmeraligner/{database}.res" %OUT_FOLDER,
    sam = "%s/{sample}/kmeraligner/{database}.sam" %OUT_FOLDER,
    seq = "%s/{sample}/kmeraligner/{database}.fsa" %OUT_FOLDER,
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

    cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} -t_db $db_path -sam 4 > {output.sam}"
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


rule Variant_identifier:
  input:
    kma_results = rules.custom_kmeralignment.output.results,
    kma_seq = rules.custom_kmeralignment.output.seq,
    indels = rules.bcftools_filter_indels.output.indels,
    indels_index = "%s/{sample}/bcftools/{database}_indels.bcf.csi" %OUT_FOLDER,
    variants = rules.bcftools_variant_call.output.variants,
    variants_index = "%s/{sample}/bcftools/{database}_variants.bcf.csi" %OUT_FOLDER
  params:
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Variant_identifier"]["options"],
    id = "{sample}"
  output:
    variant_results = "%s/{sample}/Variant_identifier/{database}.tsv" %OUT_FOLDER
  conda:
    config["analysis_settings"]["Variant_identifier"]["yaml"]
  message:
    "[Variant_identifier]: Identifying variants of {wildcards.database} on {wildcards.sample}"
  shell:
    """
    python workflow/scripts/Variant_identifier.py --sample_id {params.id}  --res {input.kma_results} --fsa {input.kma_seq} --call {input.variants} --indels {input.indels} -o {output.variant_results} {params.options}
    """