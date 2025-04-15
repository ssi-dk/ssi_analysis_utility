#------------------------------ Modules --------------------------------#

# Rule: kmeraligner
# Identifies microbial species or strain using k-mer-based alignment.
rule kmeraligner:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1]
    params:
        # Path to the kmerfinder database, KMA aligner, and taxa file.
        db_path = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kmeraligner"]["database"],
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kmeraligner"]["additional_option"],
        prefix = "{out}/{sample}/kma/{sample}"
    output:
        directory("{out}/{sample}/kma/")
    conda:
        config["analysis_settings"]["kmeraligner"]["yaml"]
    message:
        "kma  -ipe {input.R1} {input.R2} -o {params.prefix} -t_db {params.db_path} {params.add_opt}"
    shell:
        """
        # Check if the output directory exists, and skip if it does.
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1
            else
                mkdir {output}
        fi 
        # Run kmaaligner with the specified parameters and inputs.
        kma  -ipe {input.R1} {input.R2} -o {params.prefix} -t_db {params.db_path} {params.add_opt}
        """


# Rule for emm typing using BLAST
rule emm_typing:
    input:
        # Input: emm allele files and assembly file for the sample
        emm_allele_files = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["emm_typing"]["emm_allele_file"],
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]
    params:
        # Parameter for coverage threshold
        cov_per = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["emm_typing"]["cov_threshold"],
    output:
        # Output is a directory to store emm typing results for the sample
        directory("{out}/{sample}/emm_typing/")
    conda:
        config["analysis_settings"]["emm_typing"]["yaml"]
    shell:
        """
        # Check if the output directory exists, skip execution if it does
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists,skipping."
                exit 1
            else
                mkdir {output}  # Create the directory if it doesn't exist
        fi

        # Run BLAST for emm typing with specified input and output options
        blastn -query {input.emm_allele_files} \
               -subject {input.assembly} \
               -qcov_hsp_perc {params.cov_per} \
               -out {output}/blast_output.tsv \
               -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send sseq evalue bitscore'
        """

# Rule for determining assembly lineage using nucmer and delta files
rule assembly_lineage_determination:
    input:
        # Input: Reference file for lineage determination and assembly file for the sample
        reference = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["assembly_lineage_determination"]["reference_fasta_file"],
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]
    output:
        # Output is a directory to store lineage determination results
        directory("{out}/{sample}/assembly_lineage/")
    params:
        # Various parameters including delta file name, frankenfasta file, and arguments for nucmer and deltafilter
        app_path = config["analysis_settings"]["assembly_lineage_determination"]["script"],
        name = lambda wildcards: wildcards.sample,
        nucmerargs = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["assembly_lineage_determination"]["nucmerargs"],
        deltafile= lambda wildcards: wildcards.sample + ".filtered.delta",
        frankenfasta = lambda wildcards: wildcards.sample + ".frankenfasta",
        deltafilterargs= lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["assembly_lineage_determination"]["deltafilterargs"]
    conda: 
        config["analysis_settings"]["assembly_lineage_determination"]["yaml"]
    shell:
        """ 
        nucmer_path=$(which nucmer)
        deltafilder_path=$(which delta-filter)
        python {params.app_path}/convert_external_genome.py --nucmerpath $nucmer_path \
                                                            --nucmerargs {params.nucmerargs} \
                                                            --deltafilterpath $deltafilder_path \
                                                            --deltafilterargs {params.deltafilterargs} \
                                                            --reference {input.reference} \
                                                            --external {input.assembly} \
                                                            --name {params.name} \
                                                            --outputdir {output} \
        """

# Rule for Kleborate typing (pathogen typing based on assembly)
rule kleborate:
    input:
        # Input: Assembly file for the sample
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]
    params:
        # Parameter for preset options
        preset = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kleborate"]["preset"],
    output:
        # Output is a directory to store Kleborate results
        directory("{out}/{sample}/kleborate/")
    conda:
        config["analysis_settings"]["kleborate"]["yaml"]
    message:
        "kleborate -a {input.assembly}  -o {output} -p {params.preset} "
    shell:
        """
        # Check if the output directory exists, skip execution if it does
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1  # Exit without running if the directory exists
            else
                mkdir {output}  # Create the directory if it doesn't exist
        fi 
        kleborate -a {input.assembly}  -o {output} -p {params.preset}  # Run Kleborate on the assembly
        """

# Rule for CHtyper (serotype prediction based on assembly)
rule CHtyper:
    input:
        # Input: Assembly file for the sample
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]
    params:
        # Parameters for CHtyper including application path, database, threshold, coverage, and BLAST options
        app_path = config["analysis_settings"]["CHtyper"]["script"],
        database = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["CHtyper"]["database"],
        threshold = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["CHtyper"]["threshold"],
        coverage = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["CHtyper"]["coverage"],
    output:
        # Output is a directory to store CHtyper results
        directory("{out}/{sample}/chtyper/")
    conda:
        config["analysis_settings"]["CHtyper"]["yaml"]
    shell:
        """
        # Check if the output directory exists, skip execution if it does
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists, skipping."
                exit 1  # Exit without running if the directory exists
            else
                mkdir {output}  # Create the directory if it doesn't exist
        fi 
        
        blastn=$(which blastn)
        
        # Run CHtyper Python script with the specified parameters
        python {params.app_path}/CHTyper-1.0.py -i {input.assembly}  -o {output} -p {params.database} -t {params.threshold} -l {params.coverage} -b $blastn
        """
