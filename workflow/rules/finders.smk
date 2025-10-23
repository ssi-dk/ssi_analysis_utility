rule plasmidfinder:
  input:
    # Input paired-end Illumina reads.
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_PlasmidFinder.output.database
  output:
    # Output directory for plasmidfinder results.
    out_dir = directory("%s/{sample}/plasmidfinder" %output_folder),
    done = temp("%s/{sample}/plasmidfinder/plasmidfinder.done" %output_folder)
  conda:
    "../envs/plasmidfinder.yaml"
  log:
    stdout = 'Logs/{sample}/plasmidfinder.log'
  message:
    "[PlasmidFinder]: Running PlasmidFinder on {wildcards.sample}"
  shell:
    """
    mkdir -p {output.out_dir}

    cmd="plasmidfinder.py -i {input.R1} {input.R2} -o {output.out_dir} -p {input.database} -x"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
        
    echo "PlasmidFinder completed successfully for {wildcards.sample}" > {output.done}
    """

rule resfinder:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    res_database = rules.setup_ResFinder.output.database,
    point_database = rules.setup_PointFinder.output.database, #Pointfinder requires `species` definition
    disin_database = rules.setup_DisinFinder.output.database
  params:
    # Point mutation
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["resfinder"]["options"]
  output:
    out_dir = directory("%s/{sample}/resfinder" %output_folder),
    done = temp("%s/{sample}/resfinder/resfinder.done" %output_folder)
  conda:
    "../envs/resfinder.yaml"
  log:
    stdout = 'Logs/{sample}/resfinder.log'
  message:
    "[ResFinder]: Running ResFinder, PointFinder, and DisinFinder on {wildcards.sample}"
  shell:
    """
    mkdir -p {output.out_dir}

    cmd="run_resfinder.py -ifq {input.R1} {input.R2} -o {output.out_dir} --acquired -db_res {input.res_database} --disinfectant -db_disinf {input.disin_database} --point -db_point {input.point_database} {params.options}"
 
    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1

    echo "ResFinder completed successfully for {wildcards.sample}" > {output.done}
    """

rule virulencefinder:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_VirulenceFinder.output.database
  output:
    out_dir = directory("%s/{sample}/virulencefinder" %output_folder),
    done = temp("%s/{sample}/virulencefinder/virulencefinder.done" %output_folder)
  conda:
    "../envs/virulencefinder.yaml"
  log:
    stdout = 'Logs/{sample}/virulencefinder.log'
  message:
    "[VirulenceFinder]: Running VirulenceFinder on {wildcards.sample}"
  shell:
    """
    mkdir -p {output.out_dir}

    cmd="virulencefinder.py -i {input.R1} {input.R2} -o {output.out_dir} -p {input.database}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1

    echo "VirulenceFinder completed successfully for {wildcards.sample}" > {output.done}
    """

rule serotypefinder:
  input:
    R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
    R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    database = rules.setup_SerotypeFinder.output.database
  output:
    out_dir = directory("%s/{sample}/serotypefinder" %output_folder),
    done = temp("%s/{sample}/serotypefinder/serotypefinder.done" %output_folder)
  conda:
    "../envs/serotypefinder.yaml"
  log:
    stdout = 'Logs/{sample}/serotypefinder.log'
  message:
    "[SerotypeFinder]: Running SerotypeFinder on {wildcards.sample}"
  shell:
    """
    mkdir -p {output.out_dir}

    cmd="serotypefinder -i {input.R1} {input.R2} -o {output.out_dir} -p {input.database} -x"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1

    echo "SerotypeFinder completed successfully for {wildcards.sample}" > {output.done}
    """


rule spa_typing:
    input:
        assembly = rules.assembly.output.output_assembly,
        database = rules.setup_Spatyper.output.database
    output:
        done = temp("%s/{sample}/spa_typing/{assembler}.done" %output_folder)
    conda:
        "../envs/spatyper.yaml"
    log:
        stdout = 'Logs/{sample}/spatyper_{assembler}.log'
    message:
        "[Spatyping]: Running Spatyper for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        outdir=$(dirname {output.done})
        cmd="spatyper.py -i {input.assembly} -db {input.database}/ -o $outdir"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Spatyping completed successfully for {wildcards.sample} with {wildcards.assembler} assembly" > {output.done}
        """


rule amrfinder:
  input:
    assembly = rules.assembly.output.output_assembly,
    database = rules.setup_AMRFinder.output.database
  params:
    # Point mutation
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["amrfinder"]["options"]
  output:
    result = "%s/{sample}/amrfinder/{assembler}.tsv" %output_folder,
    done = temp("%s/{sample}/amrfinder/{assembler}.done" %output_folder)
  conda:
    "../envs/amrfinder.yaml"
  log:
    stdout = 'Logs/{sample}/amrfinder_{assembler}.log'
  message:
    "[AMRFinderPlus]: Running AMRFinderPlus for {wildcards.sample} using ({wildcards.assembler}) contigs"
  shell:
    """
    mkdir -p $(dirname {output.result})

    cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.options} --output {output.result}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1

    echo "AMRFinderPlus completed successfully for {wildcards.sample} with {wildcards.assembler} assembly" > {output.done}
    """


rule snp_identifier:
  input:
    kma_results = rules.custom_kmeralignment.output.results,
    variants = rules.bcftools_variant_call.output.variants,
    variants_index = rules.bcftools_variant_call.output.index,
    ref_bed = rules.fetch_genbank.output.bed,
  params:
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["snp_identifier"]["options"],
    metafile = "%s/SNP_metafile.tsv" %metadata_path
  output:
    indentified_variants = "%s/{sample}/snp_identifier/{database}.tsv" %output_folder,
    done = temp("%s/{sample}/snp_identifier/{database}.done" %output_folder)
  conda:
    "../envs/python_functions.yaml"
  log:
    stdout = "Logs/{sample}/snp_identifier_{database}.log"
  message:
    "[SNP Identifier]: Identifying SNPs of {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="python workflow/scripts/SNP_identifier.py --res {input.kma_results} --call {input.variants} --bed {input.ref_bed} --metafile {params.metafile} -o {output.indentified_variants} {params.options} > {log.stdout} 2>&1"
   
    echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1

    echo "SNP identification completed successfully of {wildcards.database} on {wildcards.sample}" > {output.done}
    """
    
rule deletion_identifier:
  input:
    kma_results = rules.custom_kmeralignment.output.results,
    kma_seq = rules.custom_kmerconsensus.output.seq,
    indels = rules.bcftools_filter_indels.output.indels,
    indels_index = rules.bcftools_filter_indels.output.index,
    variants = rules.bcftools_variant_call.output.variants,
    variants_index = rules.bcftools_variant_call.output.index,
    ref_bed = rules.fetch_genbank.output.bed,
  params:
    options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["deletion_identifier"]["options"],
    metafile = "%s/deletion_metafiles.tsv" %metadata_path
  output:
    indentified_variants = "%s/{sample}/deletion_identifier/{database}.tsv" %output_folder,
    done = temp("%s/{sample}/deletion_identifier/{database}.done" %output_folder)
  conda:
    "../envs/python_functions.yaml"
  log:
    stdout = "Logs/{sample}/Deletion_identifier_{database}.log"
  message:
    "[Deletion Identifier]: Identifying deletions of {wildcards.database} on {wildcards.sample}"
  shell:
    """
    cmd="python workflow/scripts/deletion_identifier.py --res {input.kma_results} --fsa {input.kma_seq} --call {input.variants} --indels {input.indels} --bed {input.ref_bed} --metafile {params.metafile} -o {output.indentified_variants} {params.options} > {log.stdout} 2>&1"

    echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1

    echo "Deletion identification completed successfully of {wildcards.database} on {wildcards.sample}" > {output.done}
    """

rule cdiff_repeat_identifier:
  input:
    seqs  = expand(rules.fetch_type_repeat_sequence.output.seq, TR = ["TR6", "TR10"]),
    metas = expand(rules.fetch_type_repeat_metadata.output.meta, TR = ["TR6", "TR10", "TRST"]),
    assembly = rules.assembly.output.output_assembly
  output:
    repeat_types = "%s/{sample}/CDiff_Repeat_identifier/{assembler}_repeat_types.tsv" %output_folder,
    done = temp("%s/{sample}/CDiff_Repeat_identifier/{assembler}.done" %output_folder)
  params:
    repeats = ["TR6", "TR10"],
    combos = ["TRST"]
  conda:
    "../envs/python_functions.yaml"
  log:
    stdout = "Logs/{sample}/CDiff_Repeat_identifier/{assembler}_repeat_types.log"
  message:
    "[CDiff Repeat identifier]: Identifying C. Difficile repeats in {wildcards.sample} on {wildcards.assembler} assembly"
  shell:
    """
    mkdir -p $(dirname {output.repeat_types})

    db_dir=$(dirname {input.seqs} | uniq)

    cmd="python workflow/scripts/Repeat_Identifier.py --fasta {input.assembly} --ref_seq {input.seqs} --ref_meta {input.metas} --output {output.repeat_types} --sample_id {wildcards.sample} --repeats {params.repeats} --combos {params.combos} --suffix tsv --log_file {log.stdout} > {log.stdout} 2>&1"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1 
    
    echo "C. difficile repeat identification completed successfully on {wildcards.sample} with {wildcards.assembler} assembly" > {output.done}
    """
