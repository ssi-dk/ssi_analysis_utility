rule samtools_sam_filtration:
    input:
        sam = "%s/{sample}/raw/samtools/{database}.sam" %outdir
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["view_options"]
    output:
        bam = temp("%s/{sample}/raw/samtools/{database}_filtered.bam" %outdir)
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = "%s/{sample}/custom_kmeralignment_samtools_filtration_{database}.log" %logdir
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
        bam = "%s/{sample}/raw/samtools/{database}.bam" %outdir
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["view_options"]
    output:
        bam = temp("%s/{sample}/raw/samtools/{database}_filtered.bam" %outdir)
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = "%s/{sample}/custom_kmeralignment_samtools_filtration_{database}.log" %logdir
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
        bam = "%s/{sample}/raw/samtools/{database}_filtered.bam" %outdir
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["sort_options"]
    output:
        bam_sort = temp("%s/{sample}/raw/samtools/{database}_sorted.bam" %outdir),
        index = temp("%s/{sample}/raw/samtools/{database}_sorted.bam.bai" %outdir)
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = "%s/{sample}/samtools_sort_{database}.log" %logdir
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
        reference = "%s/samtools/{database}.fasta" %database_dir
    output:
        pileup = temp("%s/{sample}/raw/bcftools/{database}_pileup.bcf" %outdir),
        index = temp("%s/{sample}/raw/bcftools/{database}_pileup.bcf.csi" %outdir)
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = "%s/{sample}/bcftools_pileup_{database}.log" %logdir
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
        options = lambda wc: sample_configs[wc.sample]["bcftools"]["view_options"]
    output:
        indels = temp("%s/{sample}/raw/bcftools/{database}_pileup_indels.bcf" %outdir),
        index = temp("%s/{sample}/raw/bcftools/{database}_pileup_indels.bcf.csi" %outdir)
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = "%s/{sample}/bcftools_filter_indels_{database}.log" %logdir
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
        variants = temp("%s/{sample}/raw/bcftools/{database}_call_variants.bcf" %outdir),
        index = temp("%s/{sample}/raw/bcftools/{database}_call_variants.bcf.csi" %outdir)
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = "%s/{sample}/bcftools_variant_call_{database}.log" %logdir
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


rule snp_identifier:
    input:
        variants = rules.bcftools_variant_call.output.variants,
        variants_index = rules.bcftools_variant_call.output.index,
    params:
        options = lambda wc: sample_configs[wc.sample]["snp_identifier"]["options"],
        metafile = "%s/SNP_metafile.tsv" %SCREENING_DIR
    output:
        indentified_variants = "%s/{sample}/raw/snp_identifier/{database}.tsv" %outdir
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/snp_identifier_{database}.log" %logdir
    message:
        "[SNP Identifier]: Identifying SNPs of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="python {SCRIPTS_DIR}/SNP_identifier.py {params.options} --call {input.variants} --metafile {params.metafile} --output {output.indentified_variants}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule deletion_identifier:
    input:
        kma_seq = rules.custom_kmerconsensus.output.seq,
        indels = rules.bcftools_filter_indels.output.indels,
        indels_index = rules.bcftools_filter_indels.output.index,
        variants = rules.bcftools_variant_call.output.variants,
        variants_index = rules.bcftools_variant_call.output.index,
        asm_aln = rules.assembly_minimap2.output.results
    params:
        options  = lambda wc: sample_configs[wc.sample]["deletion_identifier"]["options"],
        metafile = "%s/deletion_metafile.tsv" %SCREENING_DIR
    output:
        identified_variants = f"{outdir}/{{sample}}/raw/deletion_identifier/{{assembler,[^_]+}}_{{database}}.tsv" #added regex expression to ensure assemblies cannot contain '_' which our database also does
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/deletion_identifier_{assembler}_{database}.log" %logdir
    message:
        "[Deletion Identifier]: Identifying deletions of {wildcards.database} on {wildcards.sample} ({wildcards.assembler})"
    shell:
        """
        cmd="python {SCRIPTS_DIR}/deletion_identifier.py {params.options} --fsa {input.kma_seq} --call {input.variants} --mpileup {input.indels} --metafile {params.metafile} --sam {input.asm_aln} --output {output.identified_variants}"


        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule kma_filter:
    input:
        results = rules.custom_kmeralignment.output.results,
        database = rules.setup_custom_kmeraligner_index.output.names
    params:
        options = lambda wildcards: sample_configs[wildcards.sample]["kma_filter"]["options"],
        metafile = "%s/kma_filter.tsv" %SCREENING_DIR
    output:
        filtered_tsv = "%s/{sample}/raw/kma_filter/{database}.tsv" % outdir
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/KMA_results/{sample}_{database}.log" %logdir
    message:
        "[KMA Filter]: Filtering KMA .res result for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        cmd="python {SCRIPTS_DIR}/KMA_Filter.py --KMA_res {input.results} --metafile {params.metafile} --sample_id {wildcards.sample} --output {output.filtered_tsv} {params.options} > {log.stdout} 2>&1"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """


rule LREfinder:
    input:
        res = rules.custom_kmeralignment.output.results,
        matrix = rules.custom_kmeralignment.output.matrix
    # params:
    #     options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["custom_blaster"]["options"],    
    output:
        results = "%s/{sample}/raw/LREfinder/{database}.tsv" %outdir,
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/LRE-finder_{database}.log" %logdir
    message:
        "[LRE-finder]: Identify genes and mutations leading to linezolid resistance in E. faecalis and E. faecium"
    shell:
        """
        mkdir -p $(dirname {output.results})
    
        cmd="python {SCRIPTS_DIR}/LRE-Typer.py -ires {input.res} -imat {input.matrix} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "LRE-finder successfully executed" > {log.stdout}
        """


rule chtyper:
    input:
        results = rules.custom_kmeralignment.output.results
    params:
        id = 90,
        coverage = 60
    output:
        filtered_tsv = "%s/{sample}/raw/chtyper/{database}_chtyper.tsv" % outdir
    log:
        stdout = "%s/{sample}/{database}_chtyper.log" %logdir
    message:
        "[CH Typer]: Running Chtyper on {wildcards.database} assembly for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        echo "Running awk filter on {input.results}" > {log.stdout} 2>&1

        awk -F'\t' 'NR==1{{for(i=1;i<=NF;i++){{if($i=="Template_Identity")id=i;if($i=="Template_Coverage")cov=i}}print;next}} ($id+0>{params.id} && $cov+0>{params.coverage})' {input.results} > {output.filtered_tsv} 2>> {log.stdout}
        """