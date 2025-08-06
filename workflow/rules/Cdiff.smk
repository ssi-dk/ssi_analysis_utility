#------------------------------ Modules --------------------------------#
    
rule Variant_identifier:
    input:
        genotype_call_bcf = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.calls.bcf",
        indels_only_bcf = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.indels.bcf",
        genotype_call_csi = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.calls.bcf.csi",
        indels_only_csi = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.indels.bcf.csi",
    output:
        filtered_tsv = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.variants.tsv",
    params:
        input_folder = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["input_folder"],
        kma_res = lambda wildcards: os.path.join(
            output_folder,
            wildcards.sample,
            species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["input_folder"],
            "%s.res" % wildcards.sample
        ),
        kma_fsa = lambda wildcards: os.path.join(
            output_folder,
            wildcards.sample,
            species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["input_folder"],
            "%s.fsa" % wildcards.sample
        ),
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Variant_detection"]["additional_option"],
        log_dir = lambda wildcards: "%s/%s/GenotypeCalls/" % (output_folder, wildcards.sample),
        id = lambda wildcards: wildcards.sample,
        region_buffer = 5,
        overlap = 0.3,
    log:
        stdout = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.Variant_identifier.log"
    conda:
        config["analysis_settings"]["Variant_identifier"]["yaml"]
    message:
        "[Variant Identification]: Filtering KMA .res, KMA consensus .fsa, genotype calls and indels for {wildcards.sample}"
    shell:
        """
        cmd="python workflow/scripts/Variant_identifier.py \
            --sample_id {params.id} \
            --res {params.kma_res} \
            --fsa {params.kma_fsa} \
            --call {input.genotype_call_bcf} \
            --indels {input.indels_only_bcf} \
            {params.add_opt} \
            -o {output.filtered_tsv} \
            --deletion_region_buffer {params.region_buffer} \
            --partial_overlap {params.overlap}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """

rule Repeat_Identifier:
    input:
        fasta = lambda wildcards: os.path.join(
            output_folder,
            wildcards.sample,
            wildcards.assembler,
            {
                "spades": "contigs.fasta",
                "skesa": "%s.contigs.fasta" % wildcards.sample
            }[wildcards.assembler]
        ),
        repeat_fa = lambda wildcards: [
            os.path.join(
                species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["database"],
                "%s_repeat_sequences.fa" % repeats
            )
            for repeats in species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["repeats"].split()
        ]
    output:
        result = "%s/{sample}/Repeat_identifier/{assembler}_{sample}_repeat.tsv" %output_folder
    params:
        repeats = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["repeats"],
        database = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["database"],
        combo_table = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["combos"], 
        out_dir = "%s/{sample}/Repeat_identifier" %output_folder,
        sample_id = lambda wildcards: "%s" % wildcards.sample
    conda:
        config["analysis_settings"]["Repeat_identifier"]["yaml"]
    log:
        stdout = "%s/{sample}/Repeat_identifier/{assembler}_{sample}_repeat.log" %output_folder
    message:
        "[Repeat_identifier]: Running repeat analysis for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        mkdir -p {params.out_dir}

        cmd="python workflow/scripts/Repeat_Identifier.py \
            --sample_id {params.sample_id} \
            --fasta {input.fasta} \
            --repeats {params.repeats} \
            --combos {params.combo_table} \
            --db_dir {params.database} \
            --output {output.result} \
            --suffix tsv"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """
