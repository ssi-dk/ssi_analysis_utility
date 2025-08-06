#------------------------------ Modules --------------------------------#

# Rule: kmeraligner
# Identifies microbial species or strain using k-mer-based alignment.
rule kmeraligner:
    input:
        R1 = lambda wc: sample_to_illumina[wc.sample][0],
        R2 = lambda wc: sample_to_illumina[wc.sample][1],
        db_index = lambda wc: "Logs/Databases/%s.kma_index.done" % species_configs[sample_to_organism[wc.sample]]['alignment_database']['kmeraligner']['kma_index_flag'],
        db_dir = lambda wc: "%s" % getattr(rules, species_configs[sample_to_organism[wc.sample]]['alignment_database']['kmeraligner']['db']).output.database
    params:
        add_opt = lambda wc: species_configs[sample_to_organism[wc.sample]]["analyses_to_run"]["kmeraligner"]["additional_option"],
        prefix = lambda wc: "%s/%s/kmeraligner/%s" % (output_folder, wc.sample, wc.sample),
        db_prefix = lambda wc: species_configs[sample_to_organism[wc.sample]]["alignment_database"]["kmeraligner"]["kma_db_prefix"]
    output:
        res = "%s/{sample}/kmeraligner/{sample}.res" % output_folder,
        fsa = "%s/{sample}/kmeraligner/{sample}.fsa" % output_folder,
        sam = "%s/{sample}/kmeraligner/{sample}.sam" % output_folder,
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
            output_folder,
            wildcards.sample,
            species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["input_folder"],
            "%s.res" % wildcards.sample
        )
    output:
        filtered_tsv = "%s/{sample}/KMA_results/{sample}_KMA.tsv" % output_folder,
    params:
        add_opt = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["KMA_filter"]["additional_option"],
        log_dir = lambda wildcards: "%s/%s/KMA_results/" % (output_folder, wildcards.sample),
        id = lambda wildcards: "%s" % wildcards.sample
    conda:
        config["analysis_settings"]["KMA_filter"]["yaml"]
    log:
        stdout = "Logs/{sample}/KMA_results/{sample}_KMA_filter.log"
    message:
        "[KMA_filter]: Filtering KMA .res result for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        cmd="python workflow/scripts/KMAfilter.py --KMA_res {input.kma_res} --sample_id {params.id} --output {output.filtered_tsv} {params.add_opt} --log_dir {params.log_dir}"      
        
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
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
    log:
        stdout = "Logs/{sample}/emm_typing/{sample}_emm_typing.log"
    message:
        "[emm_typing]: emm_typing result for {wildcards.sample}"
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
        cmd="blastn -query {input.emm_allele_files} -subject {input.assembly} -qcov_hsp_perc {params.cov_per} -out {output}/blast_output.tsv -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send sseq evalue bitscore'"

        echo "Executing command:\n$cmd\n" >> {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
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
    log:
        stdout = "Logs/{sample}/assembly_lineage/{sample}_assembly_lineage.log"
    message:
        "[assembly_lineage_determination]: determine the assembly lineage for {wildcards.sample}"
    shell:
        """

        nucmer_path=$(which nucmer)
        deltafilder_path=$(which delta-filter)

        cmd="python {params.app_path}/convert_external_genome.py --nucmerpath $nucmer_path --nucmerargs {params.nucmerargs} --deltafilterpath $deltafilder_path --deltafilterargs {params.deltafilterargs} --reference {input.reference} --external {input.assembly} --name {params.name} --outputdir {output}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
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
    log:
        stdout = "Logs/{sample}/FilteredBAM/{sample}.{tool}.samtools_index.log"
    shell:
        """
        cmd="samtools index {input.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
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
    log:
        stdout = "Logs/{sample}/FilteredBAM/{sample}.{tool}.samtools_view.log"
    shell:
        """
        cmd="samtools view {params.add_opt} {input.sam} -o {output.filtered_bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
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
    log:
        stdout = "Logs/{sample}/FilteredBAM/{sample}.{tool}.samtools_sort.log"
    shell:
        """
        cmd="samtools sort -o {output.sorted_bam} {input.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
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
    log:
        stdout = "Logs/{sample}/GenotypeCalls/{sample}.{tool}.{tag}.bcftools_index.log"
    shell:
        """
        cmd="bcftools index -f {input.bcf}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """

rule bcftools_mpileup:
    input:
        bam = rules.samtools_sort.output.sorted_bam,
        bai = rules.samtools_index.output.bai,
        db_index = lambda wc: "Logs/Databases/%s.fa_idx.done" % species_configs[sample_to_organism[wc.sample]]['alignment_database']['kmeraligner']['fa_idx_flag'],
        db_dir = lambda wc: "%s" % getattr(rules, species_configs[sample_to_organism[wc.sample]]['alignment_database']['kmeraligner']['db']).output.database
    params:
        db_prefix = lambda wc: species_configs[sample_to_organism[wc.sample]]["alignment_database"]['kmeraligner']["kma_db_prefix"]
    output:
        mpileup = "{folder}/{sample}/GenotypeCalls/{sample}.{tool}.mpileup.bcf"
    conda:
        config["analysis_settings"]["htslib"]["yaml"]
    message:
        "[bcftools_mpileup]: Generating mpileup {input.bam} -> {output.mpileup}"
    log:
        stdout = "Logs/{sample}/GenotypeCalls/{sample}.{tool}.bcftools_mpileup.log"
    shell:
        """
        cmd="bcftools mpileup -Ob -f {input.db_dir}/{params.db_prefix}.fasta {input.bam} -o {output.mpileup}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
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
    log:
        stdout = "Logs/{sample}/GenotypeCalls/{sample}.{tool}.bcftools_call.log"
    shell:
        """
        cmd="bcftools call -mv -Ob --ploidy 1 {input.bcf} -o {output.genotypecall}"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """    

# Rule: spades_assembly
rule spades_assembly:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    output:
        contigs = "%s/{sample}/spades/contigs.fasta" %output_folder
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
        assembly = "%s/{sample}/skesa/{sample}.contigs.fasta" %output_folder
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
