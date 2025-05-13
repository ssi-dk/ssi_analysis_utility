include : "db_setups.smk"
include : "finders.smk"
include : "characterizers.smk"
include : "others.smk"


rule Legionella_summary:
    input:
        # Input: Assembly file for the sample
        legionella_sbt_folder = rules.Legionella_SBT.output,
        lag_1_blast_output_folder = rules.lag_1_detection.output
    output:
        directory("{out}/{sample}/Legionella_summary")
    conda:
        config["analysis_settings"]["SB_summarizers"]["yaml"]
    shell:
        """
        if [ ! -d {output} ];
            then
                mkdir -p {output}
        get_leg_results --legionella_sbt_file {input.legionella_sbt_folder}/legionella.sbt.tsv --lag_1_blast_output {input.lag_1_blast_output_folder}/blast.tsv --output_file {output}/Legionella.summary.tsv --sample_name {wildcards.sample}
        fi
        """



rule Spyogenes_summary:
    input:
        # Input: Assembly file for the sample
        meningotype_folder = rules.emm_typing.output
    output:
        directory("{out}/{sample}/Spyogenes_summary")
    params:
        mlst = "22"
    conda:
        config["analysis_settings"]["SB_summarizers"]["yaml"]
    shell:
        """
        if [ ! -d {output} ];
            then
                mkdir -p {output}
        get_Spyogenes_results --emm_blast_tsv {input.meningotype_folder}/blast.tsv --output_file {output}/Spyogenes.summary.tsv --sample_name {wildcards.sample}
        fi
        """
