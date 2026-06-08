rule custom_kmeralignment:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_custom_kmeraligner_index.output.names
    params:
        prefix_db = rules.setup_custom_kmeraligner_index.params.prefix    
    output:
        results = "%s/{sample}/raw/kmeraligner/{database}.res" %outdir,
        sam = temp("%s/{sample}/raw/samtools/{database}.sam" %outdir),
        matrix = temp("%s/{sample}/raw/kmeraligner/{database}.mat.gz" %outdir)
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = "%s/{sample}/custom_kmeralignment_{database}.log" %logdir
    message:
        "[kmeraligner]: Running KMA for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        SAMDIR=$(dirname {output.sam})

        mkdir -p $OUTDIR
        mkdir -p $SAMDIR

        cmd="kma -ipe {input.R1} {input.R2} -o $OUTDIR/{wildcards.database} -t_db {params.prefix_db} -na -nc -nf -sam 4 -matrix > {output.sam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

rule custom_kmerconsensus:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_custom_kmeraligner_index.output.names
    params:
        prefix_db = rules.setup_custom_kmeraligner_index.params.prefix,
    output:
        results = temp("%s/{sample}/raw/kmerconsensus/{database}.res" %outdir),
        seq = "%s/{sample}/raw/kmerconsensus/{database}.fsa" %outdir,
        aln = temp("%s/{sample}/raw/kmerconsensus/{database}.aln" %outdir)
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = "%s/{sample}/custom_kmerconsensus_{database}.log" %logdir
    message:
        "[kmerconsensus]: Running KMA for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.seq})

        
        mkdir -p $OUTDIR
        cmd="kma -ipe {input.R1} {input.R2} -o $OUTDIR/{wildcards.database} -t_db {params.prefix_db} -nf -ref_fsa"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

rule custom_bowtie2alignment:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_custom_bowtie2_index.output.bt2_1 # just locate one of the bt2 files to activate the db_setup
    params:
       options = lambda wc: sample_configs[wc.sample]["custom_bowtie2alignment"]["options"]
    output:
        sam = temp("%s/{sample}/raw/bowtie2/{database}.sam" %outdir)
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 1
    conda:
        ENVS_DIR / "bowtie2.yaml"
    log:
        stdout = "%s/{sample}/custom_bowtie2_{database}.log" %logdir
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


rule seqsero2:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        seqsero = "%s/{sample}/raw/seqsero2/SeqSero_result.tsv" %outdir
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 1
    conda:
        ENVS_DIR / "seqsero2.yaml"
    log:
        stdout = "%s/{sample}/seqsero2.log" %logdir
    message:
        "[seqsero2]: Running seqsero2 on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.seqsero})
        mkdir -p $outdir

        cmd="SeqSero2_package.py -m k -t 2 -b mem -i {input.R1} {input.R2} -d $outdir -n {wildcards.sample} -p {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


### CGE TOOLS ###
rule resfinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        res_database = rules.setup_ResFinder.output.database,
        point_database = rules.setup_PointFinder.output.database, #Pointfinder requires `species` definition
        disin_database = rules.setup_DisinFinder.output.database
    params:
        options = lambda wc: sample_configs[wc.sample]["resfinder"]["options"]
    output:
        resistance = "%s/{sample}/raw/resfinder/ResFinder_results_tab.txt" %outdir
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = "%s/{sample}/resfinder.log" %logdir
    message:
        "[ResFinder]: Running ResFinder, PointFinder, and DisinFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.resistance})
        cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o $outdir --acquired -db_res {input.res_database} --disinfectant -db_disinf {input.disin_database} --point -db_point {input.point_database} {params.options}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """



rule plasmidfinder:
    input:
        # Input paired-end Illumina reads.
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_PlasmidFinder.output.database
    output:
        replicons = "%s/{sample}/raw/plasmidfinder/results_tab.tsv" %outdir
    conda:
        ENVS_DIR / "plasmidfinder.yaml"
    log:
        stdout = "%s/{sample}/plasmidfinder.log" %logdir
    message:
        "[PlasmidFinder]: Running PlasmidFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.replicons})

        cmd="plasmidfinder.py -i {input.R1} {input.R2} -o $outdir -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1            
        """


rule virulencefinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_VirulenceFinder.output.database
    output:
        virulence = "%s/{sample}/raw/virulencefinder/results_tab.tsv" %outdir,
    conda:
        ENVS_DIR / "virulencefinder.yaml"
    log:
        stdout = "%s/{sample}/virulencefinder.log" %logdir
    message:
        "[VirulenceFinder]: Running VirulenceFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.virulence})
        cmd="python -m virulencefinder -ifq {input.R1} {input.R2} -o $outdir -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule serotypefinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_SerotypeFinder.output.database
    output:
        serotype = "%s/{sample}/raw/serotypefinder/results_tab.tsv" %outdir,
    conda:
        ENVS_DIR / "serotypefinder.yaml"
    log:
        stdout = "%s/{sample}/serotypefinder.log" %logdir
    message:
        "[SerotypeFinder]: Running SerotypeFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.serotype})
        cmd="serotypefinder -i {input.R1} {input.R2} -o $outdir -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

