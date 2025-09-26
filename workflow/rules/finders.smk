rule PlasmidFinder:
    input:
        # Input paired-end Illumina reads.
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_PlasmidFinder.output.database
    output:
        # Output directory for plasmidfinder results.
        out_dir = directory("%s/{sample}/PlasmidFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["plasmidfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/PlasmidFinder.log'
    message:
        "[PlasmidFinder]: Running PlasmidFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="plasmidfinder.py -i {input.R1} {input.R2} -o {output.out_dir} -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule ResFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        res_database = rules.setup_ResFinder.output.database,
        #point_database = rules.setup_PointFinder.output.database, #Pointfinder requires `species` definition
        disin_database = rules.setup_DisinFinder.output.database
    output:
        out_dir = directory("%s/{sample}/ResFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/ResFinder.log'
    message:
        "[ResFinder]: Running ResFinder, PointFinder, and DisinFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o {output} --acquired -db_res {input.res_database} --disinfectant -db_disinf {input.disin_database}"
 
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule VirulenceFinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_VirulenceFinder.output.database
    output:
        out_dir = directory("%s/{sample}/VirulenceFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["virulencefinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/VirulenceFinder.log'
    message:
        "[VirulenceFinder]: Running VirulenceFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="virulencefinder.py -i {input.R1} {input.R2} -o {output.out_dir} -p {input.database}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule serotypefinder:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
        database = rules.setup_SerotypeFinder.output.database
    output:
        out_dir = directory("%s/{sample}/SerotypeFinder" %OUT_FOLDER)
    conda:
        config["analysis_settings"]["serotypefinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/SerotypeFinder.log'
    message:
        "[SerotypeFinder]: Running SerotypeFinder on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="serotypefinder -i {input.R1} {input.R2} -o {output.out_dir} -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule AMRFinder:
    input:
        assembly = rules.assembly.output,
        database = rules.setup_AMRFinder.output.database
    params:
        # Point mutation
        organism = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["amrfinder"]["organism"]
    output:
        result = "%s/{sample}/AMRFinder/{assembler}.tsv" %OUT_FOLDER
    conda:
        config["analysis_settings"]["amrfinder"]["yaml"]
    log:
        stdout = 'Logs/{sample}/AMRFinder_{assembler}.log'
    message:
        "[AMRFinder]: Running AMRFinderFinder for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        mkdir -p $(dirname {output.result})

        cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.organism} --output {output.result}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule variant_identifier:
  input:
    kma_results = rules.custom_kmeralignment.output.results,
    kma_seq = rules.custom_kmerconsensus.output.seq,
    indels = rules.bcftools_filter_indels.output.indels,
    indels_index = "%s/{sample}/bcftools/{database}_indels.bcf.csi" %OUT_FOLDER,
    variants = rules.bcftools_variant_call.output.variants,
    variants_index = "%s/{sample}/bcftools/{database}_variants.bcf.csi" %OUT_FOLDER,
    ref_bed = "%s/custom/{database}.bed6" %database_path
  params:
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["Variant_identifier"]["options"]
  output:
    indentifyed_variants = "%s/{sample}/Variant_identifier/variants_{database}.tsv" %OUT_FOLDER
  conda:
    config["analysis_settings"]["Variant_identifier"]["yaml"]
  log:
    stdout = "Logs/{sample}/Variant_identifier_{database}.log"
  message:
    "[Variant Identifier]: Identifying variants of {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="python workflow/scripts/Variant_Identifier.py --sample_id {wildcards.sample}  --res {input.kma_results} --fsa {input.kma_seq} --call {input.variants} --indels {input.indels} --bed {input.ref_bed} -o {output.indentifyed_variants} {params.options}  --log_file {log.stdout} > {log.stdout} 2>&1"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """


rule CDiff_Repeat_identifier:
  input:
    seqs  = expand(rules.fetch_type_repeat_sequence.output.seq, TR = ["TR6", "TR10"]),
    metas = expand(rules.fetch_type_repeat_metadata.output.meta, TR = ["TR6", "TR10", "TRST"]),
    assembly = rules.assembly.output
  output:
    repeat_types = "%s/{sample}/CDiff_Repeat_identifier/{assembler}_repeat_types.tsv" %OUT_FOLDER
  params:
    repeats = ["TR6", "TR10"],
    combos = ["TRST"]
  conda:
    config["analysis_settings"]["Repeat_identifier"]["yaml"]
  log:
    stdout = "Logs/{sample}/CDiff_Repeat_identifier/{assembler}_repeat_types.log"
  message:
    "[CDiff Repeat identifier]: Running repeat analysis on {wildcards.assembler} assembly for {wildcards.sample}"
  shell:
    """
    mkdir -p $(dirname {output.repeat_types})

    db_dir=$(dirname {input.seqs} | uniq)

    cmd="python workflow/scripts/Repeat_Identifier.py --fasta {input.assembly} --ref_seq {input.seqs} --ref_meta {input.metas} --output {output.repeat_types} --sample_id {wildcards.sample} --repeats {params.repeats} --combos {params.combos} --suffix tsv --log_file {log.stdout} > {log.stdout} 2>&1"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1 
    """
