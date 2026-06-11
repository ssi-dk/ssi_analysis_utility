### Assembly ###

rule assembly:
    input:
        assembly = f"{outdir}/{{sample}}/raw/{{assembler}}/{{assembler}}_{{sample}}.fasta"
    output:
        assembly = f"{outdir}/{{sample}}/Assemblies/{{assembler}}_{{sample}}.fasta"
    log:
        stdout = f"{logdir}/sync_{{assembler}}_{{sample}}.log"
    message:
        "[assembly]: Syncronizing {wildcards.assembler} for {wildcards.sample} from raw location to assembly folder"
 
    shell:
        """
        mkdir -p $(dirname {output.assembly})

        cmd="rsync -av {input.assembly} {output.assembly}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 
        """

# Samtools and bcftools

rule samtools_sam_filtration:
    input:
        sam = f"{outdir}/{{sample}}/raw/samtools/{{database}}.sam"
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["view_options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/samtools/samtools_bam_filtration_{{database}}.bam")
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = f"{logdir}/samtools_sam_filtration_{{database}}_{{sample}}.log"
    message:
        "[samtools_sam_filtration]: Filtering kmeralignment output for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools view {input.sam} {params.options} -F 4 -bo {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule samtools_bam_filtration:
    input:
        bam = f"{outdir}/{{sample}}/raw/samtools/{{database}}.bam"
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["view_options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/samtools/samtools_bam_filtration_{{database}}.bam")
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = f"{logdir}/samtools_bam_filtration_{{database}}_{{sample}}.log"
    message:
        "[samtools_bam_filtration]: Filtering kmeralignment output for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools view {input.bam} {params.options} -F 4 -bo {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule samtools_sort:
    input:
        bam = f"{outdir}/{{sample}}/raw/samtools/samtools_bam_filtration_{{database}}.bam"
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["sort_options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/samtools/samtools_sort_{{database}}.bam"),
        index = temp(f"{outdir}/{{sample}}/raw/samtools/samtools_sort_{{database}}.bam.bai")
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = f"{logdir}/samtools_sort_{{database}}_{{sample}}.log"
    message:
        "[samtools_sort]: Sorting filtered bam for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools sort -o {output.results} {input.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="samtools index {output.results}"

        echo "\nIndexing Bam:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_pileup:
    input:
        bam_sort = rules.samtools_sort.output.results,
        reference = f"{database_dir}/samtools/{{database}}.fasta"
    output:
        results = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_pileup_{{database}}.bcf"),
        index = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_pileup_{{database}}.bcf.csi")
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = f"{logdir}/bcftools_pileup_{{database}}_{{sample}}.log"
    message:
        "[bcftools_pileup]: Generating mpileup for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools mpileup -Ob -f {input.reference} {input.bam_sort} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.results} -o {output.index}"

        echo "\nIndexing Pileup:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_filter_indels:
    input:
        pileup = rules.bcftools_pileup.output.results,
        pileup_index = rules.bcftools_pileup.output.index,
    params:
        options = lambda wc: sample_configs[wc.sample]["bcftools"]["view_options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_filter_indels_{{database}}.bcf"),
        index = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_filter_indels_{{database}}.bcf.csi")
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = f"{logdir}/bcftools_filter_indels_{{database}}_{{sample}}.log"
    message:
        "[bcftools_filter_indels]: Filtering indels of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools view {params.options} -Ob -o {output.results} {input.pileup}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.results} -o {output.index}"

        echo "\nIndexing Pileup:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_variant_call:
    input:
        pileup = rules.bcftools_pileup.output.results,
        pileup_index = rules.bcftools_pileup.output.index,
    output: 
        results = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_variant_call_{{database}}.bcf"),
        index = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_variant_call_{{database}}.bcf.csi")
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = f"{logdir}/bcftools_variant_call_{{database}}_{{sample}}.log"
    message:
        "[bcftools_variant_call]: Calling variant of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools call -mv -Ob --ploidy 1 {input.pileup} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.results} -o {output.index}"

        echo "\nIndexing Call:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

### SNP analysis ###

rule snp_identifier:
    input:
        variants = rules.bcftools_variant_call.output.results,
        variants_index = rules.bcftools_variant_call.output.index,
    params:
        options = lambda wc: sample_configs[wc.sample]["snp_identifier"]["options"],
        metafile = f"{SCREENING_DIR}/SNP_metafile.tsv"
    output:
        indentified_variants = f"{outdir}/{{sample}}/raw/snp_identifier/snp_identifier_{{database}}.tsv"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/snp_identifier_{{database}}_{{sample}}.log"
    message:
        "[snp_identifier]: Identifying SNPs of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="python {SCRIPTS_DIR}/SNP_identifier.py {params.options} --call {input.variants} --metafile {params.metafile} --output {output.indentified_variants}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule minimap2:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.fetch_genbank.output.fasta
    params:
        options = lambda wc: sample_configs[wc.sample]["minimap2"]["options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/minimap2/minimap2_{{assembler}}_{{database}}.sam")
    conda:
        ENVS_DIR / "minimap2.yaml"
    log:
        stdout = f"{logdir}/minimap2_{{assembler}}_{{database}}_{{sample}}.log"
    message:
        "[minimap2]: Running Minimap2 for {wildcards.database} on {wildcards.assembler} for {wildcards.sample}"
    shell:
        r"""
        mkdir -p $(dirname {output.results})

        cmd="minimap2 {params.options} {input.database} {input.assembly} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """
