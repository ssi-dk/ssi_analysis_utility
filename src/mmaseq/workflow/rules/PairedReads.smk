# Assemblers

rule spades:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = "%s/{sample}/raw/spades/{sample}.fasta" %outdir
    conda:
        ENVS_DIR / "spades.yaml"
    log:
        stdout = "%s/Assemblies/{sample}_spades.log" %logdir
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 2
    message:
        "[SPAdes]: Assemblying {wildcards.sample} using SPAdes with {threads} thread(s). This may take some time!\nInspect {log.stdout} for more details!"
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
        assembly = "%s/{sample}/raw/skesa/{sample}.fasta" %outdir
    conda:
        ENVS_DIR / "skesa.yaml"
    log:
        stdout = "%s/Assemblies/{sample}_Skesa.log" %logdir
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 2
    message:
        "[Skesa]: Assemblying {wildcards.sample} using Skesa with {threads} core(s). This may take some time!\nInspect {log.stdout} for more details!"
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
        assembly = "%s/{sample}/raw/shovill/{sample}.fasta" %outdir
    conda:
        ENVS_DIR / "shovill.yaml"
    log:
        stdout = "%s/Assemblies/{sample}_Shovill.log" %logdir
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 2
    message:
        "[Shovill]: Assemblying {wildcards.sample} using Shovill with {threads} CPU(s). This may take some time!\nInspect {log.stdout} for more details!"
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

# Custom mapping

rule kmeraligner:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_kmeraligner_index.output.names
    params:
        prefix_db = rules.setup_kmeraligner_index.params.prefix    
    output:
        results = "%s/{sample}/raw/PR/kmeraligner/{database}.res" %outdir,
        sam = temp("%s/{sample}/raw/samtools/{database}.sam" %outdir),
        matrix = temp("%s/{sample}/raw/PR/kmeraligner/{database}.mat.gz" %outdir)
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = "%s/{sample}/kmeraligner_{database}.log" %logdir
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

rule kmeraligner_consensus:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_kmeraligner_index.output.names
    params:
        prefix_db = rules.setup_kmeraligner_index.params.prefix,
    output:
        results = temp("%s/{sample}/raw/PR/kmerconsensus/{database}.res" %outdir),
        seq = "%s/{sample}/raw/PR/kmerconsensus/{database}.fsa" %outdir,
        aln = temp("%s/{sample}/raw/PR/kmerconsensus/{database}.aln" %outdir)
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = "%s/{sample}/kmeraligner_consensus_{database}.log" %logdir
    message:
        "[kmeraligner_consensus]: Running KMA for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.seq})

        
        mkdir -p $OUTDIR
        cmd="kma -ipe {input.R1} {input.R2} -o $OUTDIR/{wildcards.database} -t_db {params.prefix_db} -nf -ref_fsa"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bowtie2:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_bowtie2_index.output.bt2_1 # just locate one of the bt2 files to activate the db_setup
    params:
       options = lambda wc: sample_configs[wc.sample]["bowtie2"]["options"]
    output:
        sam = temp("%s/{sample}/raw/PR/bowtie2/{database}.sam" %outdir)
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 1
    conda:
        ENVS_DIR / "bowtie2.yaml"
    log:
        stdout = "%s/{sample}/bowtie2_{database}.log" %logdir
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

# Mappers

rule seqsero2:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        seqsero = "%s/{sample}/raw/PR/seqsero2/SeqSero_result.tsv" %outdir
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


rule resfinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        res_database = rules.setup_ResFinder.output.database
    params:
        options = lambda wc: sample_configs[wc.sample]["resfinder"]["options"]
    output:
        resistance = "%s/{sample}/raw/PR/resfinder/ResFinder_results_tab.txt" %outdir
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = "%s/{sample}/resfinder.log" %logdir
    message:
        "[resfinder]: Running ResFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.resistance})
        cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o $outdir --acquired -db_res {input.res_database} {params.options}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule pointfinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        res_database = rules.setup_ResFinder.output.database,
        point_database = rules.setup_PointFinder.output.database
    params:
        options = lambda wc: sample_configs[wc.sample]["pointfinder"]["options"]
    output:
        point_mutations = "%s/{sample}/raw/PR/pointfinder/PointFinder_results.txt" %outdir
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = "%s/{sample}/pointfinder.log" %logdir
    message:
        "[pointfinder]: Running PointFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.point_mutations})
        cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o $outdir -db_res {input.res_database} --point -db_point {input.point_database} {params.options}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule disinfinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        res_database = rules.setup_ResFinder.output.database,
        disin_database = rules.setup_DisinFinder.output.database
    params:
        options = lambda wc: sample_configs[wc.sample]["disinfinder"]["options"]
    output:
        disins = "%s/{sample}/raw/PR/disinfinder/DisinFinder_results_tab.txt" %outdir
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = "%s/{sample}/disinfinder.log" %logdir
    message:
        "[disinfinder]: Running DisinFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.disins})
        cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o $outdir -db_res {input.res_database} --disinfectant -db_disinf {input.disin_database} {params.options}"
    
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
        replicons = "%s/{sample}/raw/PR/plasmidfinder/results_tab.tsv" %outdir
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
        virulence = "%s/{sample}/raw/PR/virulencefinder/results_tab.tsv" %outdir,
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
        serotype = "%s/{sample}/raw/PR/serotypefinder/results_tab.tsv" %outdir,
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

