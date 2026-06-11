### Assembly ###

rule spades:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = f"{outdir}/{{sample}}/raw/spades/spades_{{sample}}.fasta"
    conda:
        ENVS_DIR / "spades.yaml"
    log:
        stdout = f"{logdir}/spades_{{sample}}.log"
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 2
    message:
        "[spades]: Assemblying {wildcards.sample} using SPAdes with {threads} thread(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        outdir=$(dirname {output.assembly})
        cmd="spades.py -1 {input.R1} -2 {input.R2} --threads {threads} --isolate -o $outdir"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="mv $outdir/contigs.fasta {output.assembly}"
        echo "### SPAdes Done!###\nExecuting command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule skesa:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = f"{outdir}/{{sample}}/raw/skesa/skesa_{{sample}}.fasta"
    conda:
        ENVS_DIR / "skesa.yaml"
    log:
        stdout = f"{logdir}/skesa_{{sample}}.log"
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 2
    message:
        "[skesa]: Assemblying {wildcards.sample} using Skesa with {threads} core(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        cmd="skesa --reads {input.R1},{input.R2} --contigs_out {output.assembly} --cores {threads}"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule shovill:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = f"{outdir}/{{sample}}/raw/shovill/shovill_{{sample}}.fasta"
    conda:
        ENVS_DIR / "shovill.yaml"
    log:
        stdout = f"{logdir}/shovill_{{sample}}.log"
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 2
    message:
        "[shovill]: Assemblying {wildcards.sample} using Shovill with {threads} CPU(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        mkdir -p $(dirname {output.assembly})
        outdir=$(dirname {output.assembly})

        cmd="shovill --R1 {input.R1} --R2 {input.R2} --outdir $outdir/ --force --cpus {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="mv $outdir/contigs.fa {output.assembly}"
        echo "### Shovil Done!###\nExecuting command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

### Mapping ###

rule kmeraligner:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_kmeraligner_index.output.names
    params:
        tmp_results = f"{{database}}.res",
        tmp_matrix = f"{{database}}.mat.gz",
        prefix_db = rules.setup_kmeraligner_index.params.prefix    
    output:
        results = f"{outdir}/{{sample}}/raw/PR/kmeraligner/kmeraligner_{{database}}.tsv",
        sam = temp(f"{outdir}/{{sample}}/raw/samtools/{{database}}.sam"),
        matrix = temp(f"{outdir}/{{sample}}/raw/PR/kmeraligner/kmeraligner_{{database}}.mat.gz")
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = f"{logdir}/kmeraligner_{{database}}_{{sample}}.log"
    message:
        "[kmeraligner]: Running KMA for {wildcards.database} on {wildcards.sample}"
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


rule kmeraligner_consensus:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_kmeraligner_index.output.names
    params:
        tmp_results = f"{{database}}.fsa",
        prefix_db = rules.setup_kmeraligner_index.params.prefix,
    output:
        results = temp(f"{outdir}/{{sample}}/raw/PR/kmerconsensus/kmeraligner_consensus_{{database}}.fasta")
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = f"{logdir}/kmeraligner_consensus_{{database}}_{{sample}}.log"
    message:
        "[kmeraligner_consensus]: Running KMA for {wildcards.database} on {wildcards.sample}"
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


rule bowtie2:
    input:
       R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
       R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
       database = rules.setup_bowtie2_index.output.bt2_1 # just locate one of the bt2 files to activate the db_setup
    params:
       options = lambda wc: sample_configs[wc.sample]["bowtie2"]["options"]
    output:
        sam = temp(f"{outdir}/{{sample}}/raw/PR/bowtie2/bowtie2_{{database}}.sam")
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 1
    conda:
        ENVS_DIR / "bowtie2.yaml"
    log:
        stdout = f"{logdir}/bowtie2_{{database}}_{{sample}}.log"
    message:
        "[bowtie2]: Running Bowtie2 for {wildcards.database} on {wildcards.sample} using {threads} thread(s)"
    shell:
        """
        mkdir -p $(dirname {output.sam})

        db_prefix="{input.database}"
        db_prefix="${{db_prefix%.1.bt2}}"

        cmd="bowtie2 -1 {input.R1} -2 {input.R2} -q -S {output.sam} {params.options} -x $db_prefix --threads {threads}"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

### Characterizers ###

rule seqsero2:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    params:
        tmp_results = "SeqSero_result.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/PR/seqsero2/seqsero2.tsv"
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 1
    conda:
        ENVS_DIR / "seqsero2.yaml"
    log:
        stdout = f"{logdir}/seqsero2_{{sample}}.log"
    message:
        "[seqsero2]: Running seqsero2 on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="SeqSero2_package.py -m k -t 2 -b mem -i {input.R1} {input.R2} -d $OUTDIR -n {wildcards.sample} -p {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule plasmidfinder:
    input:
        # Input paired-end Illumina reads.
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_plasmidfinder.output.database
    params:
        tmp_results = "results_tab.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/PR/plasmidfinder/plasmidfinder.tsv"
    conda:
        ENVS_DIR / "plasmidfinder.yaml"
    log:
        stdout = f"{logdir}/plasmidfinder_{{sample}}.log"
    message:
        "[plasmidfinder]: Running PlasmidFinder on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})

        cmd="plasmidfinder.py -i {input.R1} {input.R2} -o $OUTDIR -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1     

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule resfinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        res_database = rules.setup_resfinder.output.database
    params:
        tmp_results = "ResFinder_results_tab.txt",
        options = lambda wc: sample_configs[wc.sample]["resfinder"]["options"]
    output:
        results = f"{outdir}/{{sample}}/raw/PR/resfinder/resfinder.tsv"
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = f"{logdir}/resfinder_{{sample}}.log"
    message:
        "[resfinder]: Running ResFinder on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o $OUTDIR --acquired -db_res {input.res_database} {params.options}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule pointfinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        res_database = rules.setup_resfinder.output.database,
        point_database = rules.setup_pointfinder.output.database
    params:
        tmp_results = "PointFinder_results.txt",
        options = lambda wc: sample_configs[wc.sample]["pointfinder"]["options"]
    output:
        results = f"{outdir}/{{sample}}/raw/PR/pointfinder/pointfinder.tsv"
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = f"{logdir}/pointfinder_{{sample}}.log"
    message:
        "[pointfinder]: Running PointFinder on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o $OUTDIR -db_res {input.res_database} --point -db_point {input.point_database} {params.options}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule disinfinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        res_database = rules.setup_resfinder.output.database,
        disin_database = rules.setup_disinfinder.output.database
    params:
        tmp_results = "DisinFinder_results_tab.txt",
        options = lambda wc: sample_configs[wc.sample]["disinfinder"]["options"]
    output:
        results = f"{outdir}/{{sample}}/raw/PR/disinfinder/disinfinder.tsv"
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = f"{logdir}/disinfinder_{{sample}}.log"
    message:
        "[disinfinder]: Running DisinFinder on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o $OUTDIR -db_res {input.res_database} --disinfectant -db_disinf {input.disin_database} {params.options}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule virulencefinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_virulencefinder.output.database
    params:
        tmp_results = "results_tab.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/PR/virulencefinder/virulencefinder.tsv",
    conda:
        ENVS_DIR / "virulencefinder.yaml"
    log:
        stdout = f"{logdir}/virulencefinder_{{sample}}.log"
    message:
        "[virulencefinder]: Running VirulenceFinder on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        cmd="python -m virulencefinder -ifq {input.R1} {input.R2} -o $OUTDIR -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule serotypefinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_serotypefinder.output.database
    params:
        tmp_results = "results_tab.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/PR/serotypefinder/serotypefinder.tsv",
    conda:
        ENVS_DIR / "serotypefinder.yaml"
    log:
        stdout = f"{logdir}/serotypefinder_{{sample}}.log"
    message:
        "[serotypefinder]: Running SerotypeFinder on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})

        cmd="serotypefinder -i {input.R1} {input.R2} -o $OUTDIR -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """