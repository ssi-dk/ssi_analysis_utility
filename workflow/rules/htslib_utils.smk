rule samtools_sam_filtration:
    input:
        sam = "%s/{sample}/samtools/{database}.sam" %output_folder
    params:
        options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["tool_settings"]["samtools"]["view"]["options"]
    output:
        bam = temp("%s/{sample}/samtools/{database}_filtered.bam" %output_folder)
    conda:
        rule_all_functions.resolve_env(conda_envs,
                                        "htslib")
    log:
        stdout = "Logs/{sample}/custom_kmeralignment_samtools_filtration_{database}.log"
    message:
        "[custom_kmeralignment_samtools_filtration]: Filtering kmeralignment output for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools view {input.sam} {params.options} -F 4 -bo {output.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule samtools_bam_filtration:
    input:
        bam = "%s/{sample}/samtools/{database}.bam" %output_folder
    params:
        options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["tool_settings"]["samtools"]["view"]["options"]
    output:
        bam = temp("%s/{sample}/samtools/{database}_filtered.bam" %output_folder)
    conda:
        rule_all_functions.resolve_env(conda_envs,
                                        "htslib")
    log:
        stdout = "Logs/{sample}/custom_kmeralignment_samtools_filtration_{database}.log"
    message:
        "[custom_kmeralignment_samtools_filtration]: Filtering kmeralignment output for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools view {input.bam} {params.options} -F 4 -bo {output.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule samtools_sort:
    input:
        bam = "%s/{sample}/samtools/{database}_filtered.bam" %output_folder
    params:
        options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["tool_settings"]["samtools"]["sort"]["options"]
    output:
        bam_sort = temp("%s/{sample}/samtools/{database}_sorted.bam" %output_folder),
        index = temp("%s/{sample}/samtools/{database}_sorted.bam.bai" %output_folder)
    conda:
        rule_all_functions.resolve_env(conda_envs,
                                        "htslib")
    log:
        stdout = "Logs/{sample}/samtools_sort_{database}.log"
    message:
        "[samtools_sort]: Sorting filtered bam for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools sort -o {output.bam_sort} {input.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="samtools index {output.bam_sort}"

        echo "\nIndexing Bam:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_pileup:
    input:
        bam_sort = rules.samtools_sort.output.bam_sort,
        reference = "%s/samtools/{database}.fasta" %database_path
    output:
        pileup = temp("%s/{sample}/bcftools/{database}.bcf" %output_folder),
        index = temp("%s/{sample}/bcftools/{database}.bcf.csi" %output_folder)
    conda:
        rule_all_functions.resolve_env(conda_envs,
                                        "htslib")
    log:
        stdout = "Logs/{sample}/bcftools_pileup_{database}.log"
    message:
        "[bcftools_pileup]: Generating mpileup for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools mpileup -Ob -f {input.reference} {input.bam_sort} -o {output.pileup}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.pileup} -o {output.index}"

        echo "\nIndexing Pileup:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_filter_indels:
    input:
        pileup = rules.bcftools_pileup.output.pileup,
        pileup_index = rules.bcftools_pileup.output.index,
    params:
        region = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["tool_settings"]["bcftools"]["view"]["region"],
        options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["tool_settings"]["bcftools"]["view"]["options"]
    output:
        indels = temp("%s/{sample}/bcftools/{database}_indels.bcf" %output_folder),
        index = temp("%s/{sample}/bcftools/{database}_indels.bcf.csi" %output_folder)
    conda:
        rule_all_functions.resolve_env(conda_envs,
                                        "htslib")
    log:
        stdout = "Logs/{sample}/bcftools_filter_indels_{database}.log"
    message:
        "[bcftools_filter_indels]: Filtering indels of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools view -r {params.region} {params.options} -Ob -o {output.indels} {input.pileup}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.indels} -o {output.index}"

        echo "\nIndexing Pileup:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_variant_call:
    input:
        pileup = rules.bcftools_pileup.output.pileup,
        pileup_index = rules.bcftools_pileup.output.index,
    output: 
        variants = temp("%s/{sample}/bcftools/{database}_variants.bcf" %output_folder),
        index = temp("%s/{sample}/bcftools/{database}_variants.bcf.csi" %output_folder)
    conda:
        rule_all_functions.resolve_env(conda_envs,
                                        "htslib")
    log:
        stdout = "Logs/{sample}/bcftools_variant_call_{database}.log"
    message:
        "[bcftools_variant_call]: Calling variant of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools call -mv -Ob --ploidy 1 {input.pileup} -o {output.variants}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.variants} -o {output.index}"

        echo "\nIndexing Pileup:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """
