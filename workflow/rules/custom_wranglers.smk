rule kma_filter:
    input:
        results = rules.custom_kmeralignment.output.results
    output:
        filtered_tsv = "%s/{sample}/KMA_Filter/KMA_results.tsv" % OUT_FOLDER,
    params:
        options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["options"],
        log_dir = "%s/{sample}/KMA_results/" %OUT_FOLDER,
        sample = "{sample}"
    conda:
        config["analysis_settings"]["KMA_filter"]["yaml"]
    log:
        stdout = "Logs/{sample}/KMA_results/{sample}_KMA_filter.log"
    message:
        "[KMA Filter]: Filtering KMA .res result for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        python workflow/scripts/KMA_Filter.py --KMA_res {input.kma_res} --sample_id {params.sample} --output {output.filtered_tsv} {params.options} --log_dir {params.log_dir} > {log.stdout} 2>&1
        """
