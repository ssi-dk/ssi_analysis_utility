rule kma_filter:
    input:
        results = rules.kmeraligner.output.results,
        database = rules.setup_custom_kmeraligner_index.output.names
    params:
        options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kma_filter"]["options"],
        metafile = "%s/kma_filter.tsv" %metadata_path
    output:
        filtered_tsv = "%s/{sample}/kma_filter/{database}.tsv" % output_folder,
        done = temp("%s/{sample}/kma_filter/{database}.done" %output_folder)  # note folder name; pick one casing and keep it everywhere
    conda:
        "../envs/python_functions.yaml"
    log:
        stdout = "Logs/{sample}/KMA_results/{sample}_{database}.log"
    message:
        "[KMA Filter]: Filtering KMA .res result for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        cmd="python workflow/scripts/KMA_Filter.py --KMA_res {input.results} --metafile {params.metafile} --sample_id {wildcards.sample} --output {output.filtered_tsv} {params.options} > {log.stdout} 2>&1"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
    
        echo "KMA results succesfully wrangled for {wildcards.sample} on {wildcards.database}" > {output.done}
        """
