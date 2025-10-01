rule kma_filter:
    input:
        results = rules.custom_kmeralignment.output.results,
        database = rules.setup_custom_kmeraligner_index.output.names
    output:
        filtered_tsv = "%s/{sample}/KMA_Filter/KMA_results.tsv" % output_folder,
    params:
        options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["options"]
    conda:
        "../envs/python_functions.yaml"
    log:
        stdout = "Logs/{sample}/KMA_results/{sample}_{database}.log"
    message:
        "[KMA Filter]: Filtering KMA .res result for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        python workflow/scripts/KMA_Filter.py --KMA_res {input.results} --sample_id {wildcards.sample} --output {output.filtered_tsv} {params.options} --log_file {log.stdout} > {log.stdout} 2>&1
        """
