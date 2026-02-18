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
        "../envs/plasmidfinder.yaml"
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
        "../envs/resfinder.yaml"
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
        "../envs/virulencefinder.yaml"
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


rule serotypefinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_SerotypeFinder.output.database
    output:
        serotype = "%s/{sample}/serotypefinder/results_tab.tsv" %outdir,
    conda:
        "../envs/serotypefinder.yaml"
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

rule spa_typing:
    input:
        assembly = rules.assembly.output.output_assembly,
        database = rules.setup_Spatyper.output.database
    output:
        spatyper = "%s/{sample}/spatyper/{assembler}_spatype_results.tsv" %outdir
    conda:
        "../envs/spatyper.yaml"
    log:
        stdout = "%s/{sample}/spatyper_{assembler}.log" %logdir
    message:
        "[Spatyping]: Running Spatyper for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        outdir=$(dirname {output.spatyper})
        cmd="python {scripts}/SPATyper_V2.py -a {input.assembly} -d {input.database} -o {output.spatyper} -b $outdir/seq_db -l $outdir/spatyper.log "

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Spatyping completed successfully for {wildcards.sample} with {wildcards.assembler} assembly" 
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
        "../envs/amrfinder.yaml"
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
        "../envs/python_functions.yaml"
    log:
        stdout = "%s/{sample}/LRE-finder_{database}.log" %logdir
    message:
        "[LRE-finder]: Identify genes and mutations leading to linezolid resistance in E. faecalis and E. faecium"
    shell:
        """
        mkdir -p $(dirname {output.results})
    
        cmd="python {scripts}/LRE-Typer.py -ires {input.res} -imat {input.matrix} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "LRE-finder successfully executed" > {log.stdout}
        """



rule snp_identifier:
    input:
        variants = rules.bcftools_variant_call.output.variants,
        variants_index = rules.bcftools_variant_call.output.index,
    params:
        options = lambda wc: sample_configs[wc.sample]["snp_identifier"]["options"],
        metafile = "%s/SNP_metafile.tsv" %target_screening_path
    output:
        indentified_variants = "%s/{sample}/snp_identifier/{database}.tsv" %outdir
    conda:
        "../envs/python_functions.yaml"
    log:
        stdout = "%s/{sample}/snp_identifier_{database}.log" %logdir
    message:
        "[SNP Identifier]: Identifying SNPs of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="python {scripts}/SNP_identifier.py {params.options} --call {input.variants} --metafile {params.metafile} --output {output.indentified_variants}"
    
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
        metafile = "%s/deletion_metafile.tsv" %target_screening_path
    output:
        identified_variants = f"{outdir}/{{sample}}/deletion_identifier/{{assembler,[^_]+}}_{{database}}.tsv" #added regex expression to ensure assemblies cannot contain '_' which our database also does
    conda:
        "../envs/python_functions.yaml"
    log:
        stdout = "%s/{sample}/deletion_identifier_{assembler}_{database}.log" %logdir
    message:
        "[Deletion Identifier]: Identifying deletions of {wildcards.database} on {wildcards.sample} ({wildcards.assembler})"
    shell:
        """
        cmd="python {scripts}/deletion_identifier.py {params.options} --fsa {input.kma_seq} --call {input.variants} --mpileup {input.indels} --metafile {params.metafile} --sam {input.asm_aln} --output {output.identified_variants}"


        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

rule cdiff_repeat_identifier:
    input:
        seqs  = expand(rules.fetch_type_repeat_sequence.output.seq, TR = ["TR6", "TR10"]),
        metas = expand(rules.fetch_type_repeat_metadata.output.meta, TR = ["TR6", "TR10", "TRST"]),
        assembly = rules.assembly.output.output_assembly
    params:
        repeats = lambda wc: sample_configs[wc.sample]["cdiff_repeat_identifier"]["repeats"],
        combos = lambda wc: sample_configs[wc.sample]["cdiff_repeat_identifier"]["combos"]
    output:
        repeat_types = "%s/{sample}/cdiff_repeat_identifier/{assembler}_repeat_types.tsv" %outdir
    conda:
        "../envs/python_functions.yaml"
    log:
        stdout = "%s/{sample}/cdiff_repeat_identifier_{assembler}_repeat_types.log" %logdir
    message:
        "[CDiff Repeat identifier]: Identifying C. Difficile repeats in {wildcards.sample} on {wildcards.assembler} assembly"
    shell:
        """
        mkdir -p $(dirname {output.repeat_types})

        db_dir=$(dirname {input.seqs} | uniq)

        cmd="python {scripts}/Repeat_Identifier.py --fasta {input.assembly} --ref_seq {input.seqs} --ref_meta {input.metas} --output {output.repeat_types} --sample_id {wildcards.sample} --repeats {params.repeats} --combos {params.combos} --suffix tsv"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 
        """
