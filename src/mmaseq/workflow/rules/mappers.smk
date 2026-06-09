rule custom_kmeralignment:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_custom_kmeraligner_index.output.names
    params:
        tmp_results = f"{{database}}.res",
        tmp_matrix = f"{{database}}.mat.gz",
        prefix_db = rules.setup_custom_kmeraligner_index.params.prefix    
    output:
        results = f"{outdir}/{{sample}}/raw/kmeraligner/custom_kmeralignment_{{database}}.tsv",
        sam = temp(f"{outdir}/{{sample}}/raw/samtools/{{database}}.sam"),
        matrix = temp(f"{outdir}/{{sample}}/raw/kmeraligner/custom_kmeralignment_{{database}}.mat.gz")
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = f"{logdir}/custom_kmeralignment_{{database}}_{{sample}}.log"
    message:
        "[custom_kmeralignment]: Running KMA for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        SAMDIR=$(dirname {output.sam})
        mkdir -p $SAMDIR

        cmd="kma -ipe {input.R1} {input.R2} -o $OUTDIR/{wildcards.database} -t_db {params.prefix_db} -na -nc -nf -sam 4 -matrix > {output.sam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_matrix} {output.matrix} >> {log.stdout} 2>&1
        """

rule custom_kmerconsensus:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_custom_kmeraligner_index.output.names
    params:
        tmp_results = f"{{database}}.fsa",
        prefix_db = rules.setup_custom_kmeraligner_index.params.prefix,
    output:
        results = temp(f"{outdir}/{{sample}}/raw/kmerconsensus/custom_kmerconsensus_{{database}}.fasta")
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = f"{logdir}/custom_kmerconsensus_{{database}}_{{sample}}.log"
    message:
        "[custom_kmerconsensus]: Running KMA for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="kma -ipe {input.R1} {input.R2} -o $OUTDIR/{wildcards.database} -t_db {params.prefix_db} -nf -ref_fsa"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """

# rule custom_bowtie2alignment:
#     input:
#         R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
#         R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
#         database = rules.setup_custom_bowtie2_index.output.bt2_1 # just locate one of the bt2 files to activate the db_setup
#     params:
#        options = lambda wc: sample_configs[wc.sample]["custom_bowtie2alignment"]["options"]
#     output:
#         sam = temp(f"{outdir}/{{sample}}/raw/bowtie2/{{database}}.sam")
#     threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
#     priority: 1
#     conda:
#         ENVS_DIR / "bowtie2.yaml"
#     log:
#         stdout = f"{logdir}/{{sample}}/custom_bowtie2_{{database}}.log"
#     message:
#         "[custom_bowtie2alignment]: Running Bowtie2 for {wildcards.database} on {wildcards.sample} using {threads} thread(s)"
#     shell:
#         """
#         mkdir -p $(dirname {output.sam})

#         db_prefix="{input.database}"
#         db_prefix="${{db_prefix%.1.bt2}}"
        
#         cmd="bowtie2 -1 {input.R1} -2 {input.R2} -q -S {output.sam} {params.options} -x $db_prefix --threads {threads}"
#         echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
#         eval $cmd >> {log.stdout} 2>&1
#         """

        
rule assembly_minimap2:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.fetch_genbank.output.fasta
    params:
        options = lambda wc: sample_configs[wc.sample]["assembly_minimap2"]["options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/minimap2/minimap2_{{assembler}}_{{database}}.sam")
    conda:
        ENVS_DIR / "minimap2.yaml"
    log:
        stdout = f"{logdir}/minimap2_{{assembler}}_{{database}}_{{sample}}.log"
    message:
        "[assembly_minimap2]: Running Minimap2 for {wildcards.database} on {wildcards.assembler} for {wildcards.sample}"
    shell:
        r"""
        mkdir -p $(dirname {output.results})

        cmd="minimap2 {params.options} {input.database} {input.assembly} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule kma_filter:
    input:
        results = rules.custom_kmeralignment.output.results,
        database = rules.setup_custom_kmeraligner_index.output.names
    params:
        options = lambda wildcards: sample_configs[wildcards.sample]["kma_filter"]["options"],
        metafile = f"{SCREENING_DIR}/kma_filter.tsv"
    output:
        filtered_tsv = f"{outdir}/{{sample}}/raw/kma_filter/kma_filter_{{database}}.tsv"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/kma_filter_{{database}}_{{sample}}.log"
    message:
        "[KMA kma_filter]: Filtering KMA .res result for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        cmd="python {SCRIPTS_DIR}/KMA_Filter.py --KMA_res {input.results} --metafile {params.metafile} --sample_id {wildcards.sample} --output {output.filtered_tsv} {params.options} > {log.stdout} 2>&1"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """
