rule plasmidfinder:
    input:
        # Input paired-end Illumina reads.
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_PlasmidFinder.output.database
    output:
        # Output directory for plasmidfinder results.
        replicons = "%s/{sample}/plasmidfinder/results_tab.tsv" %outdir
    conda:
        ENVS_DIR / "CGE_finders.yaml"
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
        resistance = "%s/{sample}/resfinder/ResFinder_results_tab.txt" %outdir,
        tool_version = "%s/{sample}/resfinder/ResFinder_version.txt" %outdir,
    conda:
        ENVS_DIR / "CGE_finders.yaml"
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

        cmd1="run_resfinder.py --version"
        echo "Executing command:\n$cmd1\n" > {log.stdout} 2>&1
        eval $cmd1 >> {log.stdout} 2>&1

        # 2) create version file with date
        version_cmd="run_resfinder.py --version"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s\t%s\n' "ResFinder_" "$version_str" "$date_str" > {output.tool_version}
        """


rule virulencefinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_VirulenceFinder.output.database
    output:
        virulence = "%s/{sample}/virulencefinder/results_tab.tsv" %outdir,
    conda:
        ENVS_DIR / "CGE_finders.yaml"
    log:
        stdout = "%s/{sample}/virulencefinder.log" %logdir
    message:
        "[VirulenceFinder]: Running VirulenceFinder on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.virulence})
        cmd="virulencefinder.py -i {input.R1} {input.R2} -o $outdir -p {input.database} -x"

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
        result = "%s/{sample}/amrfinder/{assembler}.tsv" %outdir,
        tool_version = "%s/{sample}/amrfinder/{assembler}_amrfinder_version.txt" %outdir,
    conda:
        ENVS_DIR / "amrfinder.yaml"
    log:
        stdout = "%s/{sample}/amrfinder_{assembler}.log" %logdir
    message:
        "[AMRFinderPlus]: Running AMRFinderPlus for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        mkdir -p $(dirname {output.result})

        cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.options} --output {output.result}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # 2) create version file with date
        version_cmd="amrfinder --version"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s\t%s\n' "amrfinder_" "$version_str" "$date_str" > {output.tool_version}
        """

rule LREfinder:
    input:
        # A complete access to the wildcard is needed, if we try to call the output of different rule we have the blending of wildcards 
        res = rules.custom_kmeralignment.output.results,
        matrix = rules.custom_kmeralignment.output.matrix
    # params:
    #     options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["custom_blaster"]["options"],    
    output:
        results = "%s/{sample}/LREfinder/{database}.tsv" %outdir,
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

rule custom_blaster:
    input:
        # A complete access to the wildcard is needed, if we try to call the output of different rule we have the blending of wildcards 
        assembly = rules.assembly.output.output_assembly,
        database = rules.fetch_custom_blast_database.output.source
    params:
        options = lambda wc: sample_configs[wc.sample]["custom_blaster"]["options"]
    output:
        results = "%s/{sample}/custom_blaster/blast_{assembler}_{database}.tsv" %outdir,
        tool_version = "%s/{sample}/custom_blaster/blast_{assembler}_{database}_version.txt" %outdir,
    conda:
        ENVS_DIR / "blast.yaml"
    log:
        stdout = "%s/{sample}/custom_blaster_{assembler}_{database}.log" %logdir
    message:
        "[setup_{wildcards.database}]: Setting up the {wildcards.database} database from the temporary storage folder"
    shell:
        """
        mkdir -p $(dirname {output.results})

        cmd="blastn -subject {input.database} -query {input.assembly} -outfmt '6' -out {output.results} {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    
        # 2) create version file with date
        version_cmd=" blastn -version|head -1"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
        """

