#------------------------------ Modules --------------------------------#

# Rule: kmeraligner
# Identifies microbial species or strain using k-mer-based alignment.
rule kmeraligner:
    input:
        R1 = lambda wc: sample_to_illumina[wc.sample][0],
        R2 = lambda wc: sample_to_illumina[wc.sample][1],
        db_index = lambda wc: f"Logs/Databases/{species_configs[sample_to_organism[wc.sample]]["alignment_database"]["kmeraligner"]["kma_index_flag"]}.kma_index.done",
        db_dir = lambda wc: f"{getattr(rules, species_configs[sample_to_organism[wc.sample]]["alignment_database"]["kmeraligner"]["db"]).output.database}",
    params:
        add_opt = lambda wc: species_configs[sample_to_organism[wc.sample]]["analyses_to_run"]["kmeraligner"]["additional_option"],
        prefix = lambda wc: f"{OUT_FOLDER}/{wc.sample}/kmeraligner/{wc.sample}",
        db_prefix = lambda wc: species_configs[sample_to_organism[wc.sample]]["alignment_database"]["kmeraligner"]["kma_db_prefix"]
    output:
        res = f"{OUT_FOLDER}" + "/{sample}/kmeraligner/{sample}.res",
        fsa = f"{OUT_FOLDER}" + "/{sample}/kmeraligner/{sample}.fsa",
        sam = f"{OUT_FOLDER}" + "/{sample}/kmeraligner/{sample}.sam",
    conda:
        config["analysis_settings"]["KMA"]["yaml"]
    log:
        stdout = 'Logs/{sample}/kmeraligner.log'
    message:
        "[kmeraligner]: Running KMA on {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.res})

        cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix} {params.add_opt} -t_db {input.db_dir}/{params.db_prefix} > {output.sam}"

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

rule samtools_index:
    input:
        bam = "{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.sorted.bam"
    output:
        bai = temp("{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.sorted.bam.bai")
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[samtools_index]: Indexing {input.bam} -> {output.bai}"
    shell:
        """
        echo "Indexing {input.bam} -> {output.bai}"
        samtools index {input.bam}
        """

rule samtools_view:
    input:
        sam = "{folder}/{sample}/{tool}/{sample}.sam"
    output:
        filtered_bam = "{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.bam"
    params:
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["samtools_view"]["additional_option"]
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
        bam = "{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.bam",
    output:
        sorted_bam = "{folder}/{sample}/FilteredBAM/{sample}.{tool}.filtered.sorted.bam"
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[samtools_sort]: Sorting {input.bam} -> {output.sorted_bam}"
    shell:
        """
        echo "Sorting {input.bam} -> {output.sorted_bam}"
        samtools sort -o {output.sorted_bam} {input.bam}
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
rule bcftools_mpileup:
    input:
        bam = rules.samtools_sort.output.sorted_bam,
        bai = rules.samtools_index.output.bai,
        db_index = lambda wc: f"Logs/Databases/{species_configs[sample_to_organism[wc.sample]]['alignment_database']['fa_idx_flag']}.fa_idx.done",
        db_dir = lambda wc: f"{getattr(rules, species_configs[sample_to_organism[wc.sample]]['alignment_database']['db']).output.database}",
    params:
        db_prefix = lambda wc: species_configs[sample_to_organism[wc.sample]]["alignment_database"]["kma_db_prefix"]
    output:
        mpileup = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf"
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[bcftools_mpileup]: Generating mpileup {input.bam} -> {output.mpileup}"
    shell:
        """
        bcftools mpileup -Ob -f {input.db_dir}/{params.db_prefix}.fasta {input.bam} -o {output.mpileup}
        """
        
# Rule: bcftools_view_filter
rule bcftools_view_filter:
    input:
        bcf = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf",
        csi = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf.csi",  # ensure indexing happens
    output:
        indels_only = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.indels.bcf",
    params:
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

        cmd="bcftools view -r {params.region} {params.add_opt} -Ob -o {output.indels_only} {input.bcf}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

rule bcftools_call:
    input:
        bcf = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf",
        csi = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf.csi"
    output:
        genotypecall = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.calls.bcf"
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[bcftools_call]: Calling genotypes from {input.bcf} -> {output.genotypecall}"
    shell:
        """
        bcftools call -mv -Ob --ploidy 1 {input.bcf} -o {output.genotypecall}
        """    
    
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
        python workflow/scripts/Variant_identifier.py --sample_id {params.id}  --res {params.kma_res} --fsa {params.kma_fsa} --call {input.genotype_call_bcf} --indels {input.indels_only_bcf} {params.add_opt} -o {output.filtered_tsv} --deletion_region_buffer {params.region_buffer} --partial_overlap {params.overlap}
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
        ),
        repeat_fa = lambda wildcards: [
            os.path.join(
                species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["database"],
                f"{repeats}_repeat_sequences.fa"
            )
            for repeats in species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Repeat_identifier"]["repeats"].split()
        ]
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

rule Salmonella_seqsero:
    input:
        R1 = lambda wc: sample_to_illumina[wc.sample][0],
        R2 = lambda wc: sample_to_illumina[wc.sample][1],
        flag = "Logs/Databases/seqsero2_db.done"
    output:
        txt = "%s/{sample}/seqsero2/SeqSero_result.txt" %OUT_FOLDER,
        tsv = "%s/{sample}/seqsero2/SeqSero_result.tsv" %OUT_FOLDER
    threads:
        min(workflow.cores, 8)
    conda:
        config["analysis_settings"]["seqsero2"]["yaml"]
    params:
        sample_id = lambda wildcards: f"{wildcards.sample}",
        out_dir = "%s/{sample}/seqsero2" %OUT_FOLDER,
    log:
        stdout = "%s/{sample}/seqsero2/SeqSero_log.txt" %OUT_FOLDER
    message:
        "[Salmonella_seqsero]: Predict Salmonella serovar"
    shell:
        """
        mkdir -p $(dirname {output.txt})

        SeqSero2_package.py -m k -t 2 -i {input.R1} {input.R2} -d {params.out_dir} -n {params.sample_id} -p {threads} -b mem > {log.stdout} 2>&1
        """
