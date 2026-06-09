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
        "[plasmidfinder]: Running PlasmidFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.replicons})

        cmd="plasmidfinder.py -i {input.R1} {input.R2} -o $outdir -p {input.database} -x"

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
        resistance = "%s/{sample}/raw/resfinder/ResFinder_results_tab.txt" %outdir
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
        point_mutations = "%s/{sample}/raw/pointfinder/PointFinder_results.txt" %outdir
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
        disins = "%s/{sample}/raw/disinfinder/DisinFinder_results_tab.txt" %outdir
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
        "[virulencefinder]: Running VirulenceFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.virulence})
        cmd="python -m virulencefinder -ifq {input.R1} {input.R2} -o $outdir -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule amrfinder:
    input:
        assembly = rules.assembly.output.output_assembly,
        database = rules.setup_AMRFinder.output.database
    params:
        options = lambda wc: sample_configs[wc.sample]["amrfinder"]["options"]
    output:
        result = "%s/{sample}/raw/amrfinder/{assembler}.tsv" %outdir
    conda:
        ENVS_DIR / "amrfinder.yaml"
    log:
        stdout = "%s/{sample}/amrfinder_{assembler}.log" %logdir
    message:
        "[amrfinder]: Running AMRFinderPlus for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        mkdir -p $(dirname {output.result})

        cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.options} --output {output.result}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
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
        "[LRE-finder]: Identify genes and mutations for linezolid resistance in {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.results})
    
        cmd="python {SCRIPTS_DIR}/LRE-Typer.py -ires {input.res} -imat {input.matrix} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "LRE-finder successfully executed" > {log.stdout}
        """

rule custom_blaster:
    input:
        # A complete access to the wildcard is needed, if we try to call the output of different rule we have the blending of wildcards 
        assembly = rules.assembly.output.output_assembly,
        database = rules.fetch_custom_blast_database.output.source
    params:
        options = lambda wc: sample_configs[wc.sample]["custom_blaster"]["options"]
    output:
        results = "%s/{sample}/raw/custom_blaster/blast_{assembler}_{database}.tsv" %outdir
    conda:
        ENVS_DIR / "blast.yaml"
    log:
        stdout = "%s/{sample}/custom_blaster_{assembler}_{database}.log" %logdir
    message:
        "[custom_blaster]: Blasting {wildcards.database} against {wildcards.sample} database from the temporary storage folder"
    shell:
        """
        mkdir -p $(dirname {output.results})

        cmd="blastn -subject {input.database} -query {input.assembly} -outfmt '6' -out {output.results} {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

