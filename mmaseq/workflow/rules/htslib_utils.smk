rule samtools_sam_filtration:
    input:
        sam = "%s/{sample}/samtools/{database}.sam" %OUTDIR
    params:
        options = lambda wc: SAMPLE_CONFIGS[wc.sample]["samtools"]["view_options"]
    output:
        bam = temp("%s/{sample}/samtools/{database}_filtered.bam" %OUTDIR)
    conda:
        "../envs/htslib.yaml"
    log:
        stdout = "%s/{sample}/custom_kmeralignment_samtools_filtration_{database}.log" %LOGDIR
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
        bam = "%s/{sample}/samtools/{database}.bam" %OUTDIR
    params:
        options = lambda wc: SAMPLE_CONFIGS[wc.sample]["samtools"]["view_options"]
    output:
        bam = temp("%s/{sample}/samtools/{database}_filtered.bam" %OUTDIR)
    conda:
        "../envs/htslib.yaml"
    log:
        stdout = "%s/{sample}/custom_kmeralignment_samtools_filtration_{database}.log" %LOGDIR
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
        bam = "%s/{sample}/samtools/{database}_filtered.bam" %OUTDIR
    params:
        options = lambda wc: SAMPLE_CONFIGS[wc.sample]["samtools"]["sort_options"]
    output:
        bam_sort = temp("%s/{sample}/samtools/{database}_sorted.bam" %OUTDIR),
        index = temp("%s/{sample}/samtools/{database}_sorted.bam.bai" %OUTDIR)
    conda:
        "../envs/htslib.yaml"
    log:
        stdout = "%s/{sample}/samtools_sort_{database}.log" %LOGDIR
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
        reference = "%s/samtools/{database}.fasta" %DATABASE_PATH
    output:
        pileup = temp("%s/{sample}/bcftools/{database}_pileup.bcf" %OUTDIR),
        index = temp("%s/{sample}/bcftools/{database}_pileup.bcf.csi" %OUTDIR)
    conda:
        "../envs/htslib.yaml"
    log:
        stdout = "%s/{sample}/bcftools_pileup_{database}.log" %LOGDIR
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
        options = lambda wc: SAMPLE_CONFIGS[wc.sample]["bcftools"]["view_options"]
    output:
        indels = temp("%s/{sample}/bcftools/{database}_pileup_indels.bcf" %OUTDIR),
        index = temp("%s/{sample}/bcftools/{database}_pileup_indels.bcf.csi" %OUTDIR)
    conda:
        "../envs/htslib.yaml"
    log:
        stdout = "%s/{sample}/bcftools_filter_indels_{database}.log" %LOGDIR
    message:
        "[bcftools_filter_indels]: Filtering indels of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools view {params.options} -Ob -o {output.indels} {input.pileup}"

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
        variants = temp("%s/{sample}/bcftools/{database}_call_variants.bcf" %OUTDIR),
        index = temp("%s/{sample}/bcftools/{database}_call_variants.bcf.csi" %OUTDIR)
    conda:
        "../envs/htslib.yaml"
    log:
        stdout = "%s/{sample}/bcftools_variant_call_{database}.log" %LOGDIR
    message:
        "[bcftools_variant_call]: Calling variant of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools call -mv -Ob --ploidy 1 {input.pileup} -o {output.variants}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.variants} -o {output.index}"

        echo "\nIndexing Call:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """