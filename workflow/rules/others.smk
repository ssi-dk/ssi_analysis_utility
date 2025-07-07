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
        "[kmeraligner]: Running KMA on {wildcards.sample}"
    shell:
        """
        mkdir -p {output}
        
        cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} {params.add_opt} -t_db {input.database}/{params.db_prefix}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# Rule: kmeraligner
# Identifies microbial species or strain using k-mer-based alignment.
rule Cdiff_KMA_Toxin:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_CdiffToxin.output.database
    params:
        db_prefix = "Cdiff_Toxin",
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Cdiff_KMA_Toxin"]["additional_option"],
        prefix = "%s/{sample}/Cdiff_KMA_Toxin/{sample}" % OUT_FOLDER
    output:
        aln = temp("%s/{sample}/Cdiff_KMA_Toxin/{sample}.aln" % OUT_FOLDER),
        res = "%s/{sample}/Cdiff_KMA_Toxin/{sample}.res" % OUT_FOLDER,
        fsa = "%s/{sample}/Cdiff_KMA_Toxin/{sample}.fsa" % OUT_FOLDER,
        sam = temp("%s/{sample}/Cdiff_KMA_Toxin/{sample}.sam" % OUT_FOLDER),
        vcf_gz = temp("%s/{sample}/Cdiff_KMA_Toxin/{sample}.vcf.gz" % OUT_FOLDER),
    conda:
        config["analysis_settings"]["Clostridioides_difficile_db"]["yaml"]
    log:
        stdout = 'Logs/{sample}/Cdiff_KMA_Toxin.log'
    message:
        "[Cdiff_KMA_Toxin]: Running Clostridium difficile kmer aligner on {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.sam})

        cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} {params.add_opt} -t_db {input.database}/{params.db_prefix} > {output.sam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

rule KMA_filter:
    input:
        kma_res = lambda wildcards: os.path.join(
            OUT_FOLDER,
            wildcards.sample,
            species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["input_folder"],
            f"{wildcards.sample}.res"
        )
    output:
        filtered_tsv = f"{OUT_FOLDER}" + "/{sample}/KMA_results/{sample}_KMA.tsv"
    params:
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["additional_option"],
        log_dir = lambda wildcards: f"{OUT_FOLDER}/{wildcards.sample}/KMA_results/",
        id = lambda wildcards: f"{wildcards.sample}",
    conda:
        config["analysis_settings"]["KMA_filter"]["yaml"]
    log:
        stdout = "Logs/{sample}/KMA_results/{sample}_KMA_filter.log"
    message:
        "[KMA_filter]: Filtering KMA .res result for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        python workflow/scripts/KMAfilter.py --KMA_res {input.kma_res} --sample_id {params.id} --output {output.filtered_tsv} {params.add_opt} --log_dir {params.log_dir} > {log.stdout} 2>&1
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

# # Rule for Kleborate typing (pathogen typing based on assembly)
# rule kleborate:
#     input:
#         # Input: Assembly file for the sample
#         assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]
#     params:
#         # Parameter for preset options
#         preset = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kleborate"]["preset"],
#     output:
#         # Output is a directory to store Kleborate results
#         directory("{out}/{sample}/kleborate/")
#     conda:
#         config["analysis_settings"]["kleborate"]["yaml"]
#     message:
#         "kleborate -a {input.assembly}  -o {output} -p {params.preset} "
#     shell:
#         """
#         # Check if the output directory exists, skip execution if it does
#         if [ -d {output} ]; 
#             then
#                 echo "Directory {output} exists, skipping."
#                 exit 1  # Exit without running if the directory exists
#             else
#                 mkdir {output}  # Create the directory if it doesn't exist
#         fi 
#         kleborate -a {input.assembly}  -o {output} -p {params.preset}  # Run Kleborate on the assembly
#         """

# # Rule for CHtyper (serotype prediction based on assembly)
# rule CHtyper:
#     input:
#         # Input: Assembly file for the sample
#         assembly = lambda wildcards: sample_to_assembly_file[wildcards.sample]
#     params:
#         # Parameters for CHtyper including application path, database, threshold, coverage, and BLAST options
#         app_path = config["analysis_settings"]["CHtyper"]["script"],
#         database = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["CHtyper"]["database"],
#         threshold = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["CHtyper"]["threshold"],
#         coverage = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["CHtyper"]["coverage"],
#     output:
#         # Output is a directory to store CHtyper results
#         directory("{out}/{sample}/chtyper/")
#     conda:
#         config["analysis_settings"]["CHtyper"]["yaml"]
#     shell:
#         """
#         # Check if the output directory exists, skip execution if it does
#         if [ -d {output} ]; 
#             then
#                 echo "Directory {output} exists, skipping."
#                 exit 1  # Exit without running if the directory exists
#             else
#                 mkdir {output}  # Create the directory if it doesn't exist
#         fi 
        
#         blastn=$(which blastn)
        
#         # Run CHtyper Python script with the specified parameters
#         python {params.app_path}/CHTyper-1.0.py -i {input.assembly}  -o {output} -p {params.database} -t {params.threshold} -l {params.coverage} -b $blastn
#         """

# Rule: samtools filter
rule samtools_view:
    input:
        sam = "{folder}/{sample}/{tool}/{sample}.sam"
    output:
        filtered_bam = temp("{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.bam")
    params:
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["samtools_view"]["additional_option"],
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[samtools_view]: Filtering {input.sam} to {output.filtered_bam}"
    shell:
        """
        echo "Filtering {input.sam} -> {output.filtered_bam}"
        samtools view {params.add_opt} {input.sam} -o {output.filtered_bam}
        """

rule samtools_sort:
    input:
        filtered_bam = "{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.bam"
    output:
        sorted_bam = temp("{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.sorted.bam")
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[samtools_sort]: Sorting {input.filtered_bam} -> {output.sorted_bam}"
    shell:
        """
        echo "Sorting {input.filtered_bam} -> {output.sorted_bam}"
        samtools sort -o {output.sorted_bam} {input.filtered_bam}
        """

rule samtools_index:
    input:
        sorted_bam = "{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.sorted.bam"
    output:
        bam_bai = temp("{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.sorted.bam.bai")
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[samtools_index]: Indexing {input.sorted_bam} -> {output.bam_bai}"
    shell:
        """
        echo "Indexing {input.sorted_bam} -> {output.bam_bai}"
        samtools index {input.sorted_bam}
        """

# Rule: bcftools_genotypecall
# Identifies microbial species or strain using k-mer-based alignment.

rule bcftools_mpileup:
    input:
        sorted_bam = "{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.sorted.bam",
        database = rules.setup_CdiffToxin.output.database
    params:
        db_prefix = rules.setup_CdiffToxin.params.db_toxin
    output:
        mpileup = temp("{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf")
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[bcftools_mpileup]: Generating mpileup {input.sorted_bam} -> {output.mpileup}"
    shell:
        """
        bcftools mpileup -Ob -f {input.database}/{params.db_prefix}.fasta {input.sorted_bam} -o {output.mpileup}
        """
        
# Rule: bcftools_view_filter
rule bcftools_view_filter:
    input:
        bcf = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf",
        database = rules.setup_CdiffToxin.output.database,
    output:
        indels_only = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.indels.bcf",
        bcf_idx = temp("{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf.csi")
    params:
        db_prefix = rules.setup_CdiffToxin.params.db_toxin,
        region = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["bcftools_view_filter"]["region"],
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["bcftools_view_filter"]["additional_option"],
    conda:
        config["analysis_settings"]["htslib"]["yaml"],
    log:
        stdout = '{folder}/{sample}/GenotypeCalls/{sample}.{tool}.bcftools_filter.log'
    message:
        "[bcftools_view_filter]: Filtering the mpileup from {input.bcf}",
    shell:
        """
        # Check if the index exists, and if not, index the BCF file first
        if [ ! -f {input.bcf}.csi ]; then
            echo "Indexing {input.bcf}"
            bcftools index {input.bcf} -o {output.bcf_idx}
        fi

        cmd="bcftools view -r {params.region} {params.add_opt} -Ob -o {output.indels_only} {input.bcf}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_call:
    input:
        mpileup = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf"
    output:
        genotypecall = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.calls.bcf"
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[bcftools_call]: Calling genotypes from {input.mpileup} -> {output.genotypecall}"
    shell:
        """
        bcftools call -mv -Ob --ploidy 1 {input.mpileup} -o {output.genotypecall}
        """

rule bcftools_index:
    input:
        bcf = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.{tag}.bcf"
    output:
        csi = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.{tag}.bcf.csi"
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[bcftools_index]: Indexing {input.bcf} -> {output.csi}"
    shell:
        """
        echo "Indexing {input.bcf} -> {output.csi}"
        bcftools index -f {input.bcf}
        """

rule Variant_identifier:
    input:
        genotype_call = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.calls.bcf",
        indels_only = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.indels.bcf"
    output:
        filtered_tsv = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.variants.tsv"
    params:
        input_folder = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["input_folder"],
        kma_res = lambda wildcards: os.path.join(
            OUT_FOLDER,
            wildcards.sample,
            species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["input_folder"],
            f"{wildcards.sample}.res"
        ),
        kma_fsa = lambda wildcards: os.path.join(
            OUT_FOLDER,
            wildcards.sample,
            species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["input_folder"],
            f"{wildcards.sample}.fsa"
        ),
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Variant_detection"]["additional_option"],
        log_dir = lambda wildcards: f"{OUT_FOLDER}/{wildcards.sample}/GenotypeCalls/",
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
        python workflow/scripts/Cdiff_wrangler_variant.py --sample_id {params.id}  --res {params.kma_res} --fsa {params.kma_fsa} --call {input.genotype_call} --indels {input.indels_only} {params.add_opt} -o {output.filtered_tsv} --deletion_region_buffer {params.region_buffer} --partial_overlap {params.overlap}
        """

# Rule: spades_assembly
rule spades_assembly:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    output:
        contigs = "%s/{sample}/spades/contigs.fasta" %OUT_FOLDER
    conda:
        config["analysis_settings"]["spades"]["yaml"]
    log:
        stdout = "Logs/{sample}/spades.log"
    message:
        "[spades_assembly]: Perform assembly using spades on {wildcards.sample}, this will take some time!"
    threads:
        min(workflow.cores, 8)
    shell:
        """
        cmd="spades.py -1 {input.R1} -2 {input.R2} --threads {threads} --isolate --only-assembler -o $(dirname {output.contigs})"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# Rule: Skesa_assembly
#  DeBruijn graph-based de-novo assembler for microbial genomes
rule skesa_assembly:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    output:
        assembly = "%s/{sample}/skesa/{sample}.contigs.fasta" %OUT_FOLDER
    conda:
        config["analysis_settings"]["skesa"]["yaml"]
    log:
        stdout = "Logs/{sample}/skesa.log"
    message:
        "[skesa_assembly]: Perform assembly using skesa on {wildcards.sample}, this will take some time!"
    threads:
        min(workflow.cores, 8)
    shell:
        """
        mkdir -p $(dirname {output.assembly})
        
        cmd="skesa --reads {input.R1},{input.R2} --contigs_out {output.assembly} --cores {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """      

rule Repeat_Identifier:
    input:
        fasta = lambda wildcards: os.path.join(
            OUT_FOLDER,
            wildcards.sample,
            wildcards.assembler,
            {
                "spades": "contigs.fasta",
                "skesa": f"{wildcards.sample}.contigs.fasta"
            }[wildcards.assembler]
        )
    output:
        result = "%s/{sample}/Repeat_identifier/{assembler}_{sample}_repeat.tsv" %OUT_FOLDER
    params:
        repeats = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["repeats"],
        database = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["database"],
        combo_table = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["combos"], 
        out_dir = "%s/{sample}/Repeat_identifier" %OUT_FOLDER,
        sample_id = lambda wildcards: f"{wildcards.sample}",
    conda:
        config["analysis_settings"]["Repeat_identifier"]["yaml"]
    log:
        stdout = "%s/{sample}/Repeat_identifier/{assembler}_{sample}_repeat.log" %OUT_FOLDER
    message:
        "[Repeat_identifier]: Running repeat analysis for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        mkdir -p {params.out_dir}

        python workflow/scripts/Repeat_Identifier.py --sample_id {params.sample_id} --fasta {input.fasta} --repeats {params.repeats} --combos {params.combo_table} --db_dir {params.database} --output {output.result} --suffix tsv > {log.stdout} 2>&1
        """
