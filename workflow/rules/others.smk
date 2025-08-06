#------------------------------ Modules --------------------------------#

# Rule: kmeraligner
# Identifies microbial species or strain using k-mer-based alignment.
rule kmeraligner:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_EcoliKmerAligner.output.database
    params:
        # Path to the kmerfinder database, KMA aligner, and taxa file.
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kmeraligner"]["additional_option"],
        db_prefix = rules.setup_EcoliKmerAligner.params.db_prefix,
        prefix = "%s/{sample}/kmeraligner/{sample}" %OUT_FOLDER
    output:
        directory("%s/{sample}/kmeraligner/" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["kmeraligner"]["yaml"]
    log:
        stdout = 'Logs/{sample}/kmeraligner.log'
    message:
        "[EcoliKmerAligner]: Running EcoliKmerAligner on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}
        
        cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} {params.add_opt} -t_db {input.database}/{params.db_prefix}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """
