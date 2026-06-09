rule plasmidfinder:
    input:
        # Input paired-end Illumina reads.
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_plasmidfinder.output.database
    params:
        tmp_results = "results_tab.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/plasmidfinder/plasmidfinder.tsv"
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
        results = f"{outdir}/{{sample}}/raw/resfinder/resfinder.tsv"
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
        results = f"{outdir}/{{sample}}/raw/pointfinder/pointfinder.tsv"
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
        results = f"{outdir}/{{sample}}/raw/disinfinder/disinfinder.tsv"
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
        results = f"{outdir}/{{sample}}/raw/virulencefinder/virulencefinder.tsv",
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


rule amrfinder:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.setup_amrfinder.output.database
    params:
        options = lambda wc: sample_configs[wc.sample]["amrfinder"]["options"]
    output:
        result = f"{outdir}/{{sample}}/raw/amrfinder/amrfinder_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "amrfinder.yaml"
    log:
        stdout = f"{logdir}/amrfinder_{{assembler}}_{{sample}}.log"
    message:
        "[amrfinder]: Running AMRFinderPlus for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR
        
        cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.options} --output {output.result}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

rule lrefinder:
    input:
        res = rules.custom_kmeralignment.output.results,
        matrix = rules.custom_kmeralignment.output.matrix
    # params:
    #     options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["custom_blaster"]["options"],    
    output:
        results = f"{outdir}/{{sample}}/raw/lrefinder/lrefinder_{{database}}.tsv",
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/lrefinder_{{database}}_{{sample}}.log"
    message:
        "[LREfinder]: Identify genes and mutations for linezolid resistance in {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="python {SCRIPTS_DIR}/LRE-Typer.py -ires {input.res} -imat {input.matrix} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

rule custom_blaster:
    input:
        # A complete access to the wildcard is needed, if we try to call the output of different rule we have the blending of wildcards 
        assembly = rules.assembly.output.assembly,
        database = rules.fetch_custom_blast_database.output.source
    params:
        options = lambda wc: sample_configs[wc.sample]["custom_blaster"]["options"]
    output:
        results = f"{outdir}/{{sample}}/raw/custom_blaster/custom_blaster_{{assembler}}_{{database}}.tsv"
    conda:
        ENVS_DIR / "blast.yaml"
    log:
        stdout = f"{logdir}/custom_blaster_{{assembler}}_{{database}}_{{sample}}.log"
    message:
        "[custom_blaster]: Blasting {wildcards.database} against {wildcards.sample} database from the temporary storage folder"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="blastn -subject {input.database} -query {input.assembly} -outfmt '6' -out {output.results} {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

