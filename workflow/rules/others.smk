#------------------------------ Modules --------------------------------#

# Rule for resistance gene detection using BLAST
rule resistance_gene_detection:
    input:
        # Input: AMR gene database file and assembly file for the sample
        amr_genes = config["analysis_settings"]["resistance_gene_detection"]["query_fasta_resistance_gene_detection"],
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]  # Use wildcards to get the sample's assembly file
    params:
        # Parameters for BLAST, including identity and coverage thresholds
        id_per = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["resistance_gene_detection"]["pident_threshold"],
        cov_per = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["resistance_gene_detection"]["cov_threshold"],
    output:
        # Output is a directory to store BLAST results for the sample
        directory("{out}/{sample}/blast/")
    shell:
        """
        # Check if the output directory exists, skip execution if it does
        if [ -d {output} ]; 
            then
                echo "Directory {output} exists,skipping."
                exit 1  # Exit without running if the directory exists
            else
                mkdir {output}  # Create the directory if it doesn't exist
        fi

        # Run BLAST with specified input and output options
        blastn -query {input.amr_genes} \
               -subject {input.assembly} \
               -out {output}/blast_output.tsv \
               -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send sseq evalue bitscore' \
               -perc_identity {params.id_per} \
               -qcov_hsp_perc {params.cov_per}
        """

# Rule for emm typing using BLAST
rule emm_typing:
    input:
        # Input: emm allele files and assembly file for the sample
        emm_allele_files = config["analysis_settings"]["emm_typing"]["emm_allele_file"],
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]
    params:
        # Parameter for coverage threshold
        cov_per = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["emm_typing"]["cov_threshold"],
    output:
        # Output is a directory to store emm typing results for the sample
        directory("{out}/{sample}/emm_typing/")
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
        reference = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["assembly_lineage_determination"]["reference_fasta_file"],
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]
    output:
        # Output is a directory to store lineage determination results
        directory("{out}/{sample}/assembly_lineage/")
    params:
        # Various parameters including delta file name, frankenfasta file, and arguments for nucmer and deltafilter
        name = lambda wildcards: wildcards.sample,
        deltafile= lambda wildcards: wildcards.sample + ".filtered.delta",
        frankenfasta = lambda wildcards: wildcards.sample + ".frankenfasta",
        nucmerargs=lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["assembly_lineage_determination"]["nucmerargs"],
        nucmerpath=config["analysis_settings"]["assembly_lineage_determination"]["nucmerpath"],
        deltafilterpath=config["analysis_settings"]["assembly_lineage_determination"]["deltafilterpath"],
        deltafilterargs= lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["assembly_lineage_determination"]["deltafilterargs"]
    run:
        # Create output directory if it doesn't exist, and run the lineage determination workflow
        if not os.path.exists(str(output)):
            os.mkdir(str(output))
            external_genome = convert_external_genome.Genome()
            external_genome.import_fasta_file(str(input.assembly))  # Import assembly as external genome
            convert_external_genome.generate_delta_file(params.nucmerpath, 
                                                        params.nucmerargs, 
                                                        params.deltafilterpath,
                                                        params.deltafilterargs, 
                                                        params.name, 
                                                        input.reference,
                                                        input.assembly, 
                                                        str(output) )  # Generate delta file using nucmer
            franken_genome = convert_external_genome.Genome()
            convert_external_genome.parse_delta_file((os.path.join(str(output), params.deltafile)),
                                                    franken_genome, 
                                                    external_genome)  # Parse delta file to create Frankenfasta
            franken_genome.write_to_fasta_file( (os.path.join(str(output), str(params.frankenfasta))), str(params.name) + " ref:")  # Write the Frankenfasta file
        else:
            print("Directory {output} exists,skipping.")
            exit()  # Skip if directory already exists


# Rule for Kleborate typing (pathogen typing based on assembly)
rule kleborate:
    input:
        # Input: Assembly file for the sample
        assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]
    params:
        # Parameter for preset options
        preset = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["kleborate"]["preset"],
    output:
        # Output is a directory to store Kleborate results
        directory("{out}/{sample}/kleborate/")
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
        app_path = config["analysis_settings"]["chtyper"]["path"],
        database = config["analysis_settings"]["chtyper"]["database"],
        threshold = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["chtyper"]["threshold"],
        coverage = lambda wildcards: config["Species"][sample_to_organism[wildcards.sample]]["analyses_to_run"]["chtyper"]["coverage"],
        blast = config["analysis_settings"]["chtyper"]["blast"]
    output:
        # Output is a directory to store CHtyper results
        directory("{out}/{sample}/chtyper/")
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
        # Run CHtyper Python script with the specified parameters
        python {params.app_path}/CHTyper-1.0.py -i {input.assembly}  -o {output} -p {params.database} -t {params.threshold} -l {params.coverage} -b {params.blast}
        """
