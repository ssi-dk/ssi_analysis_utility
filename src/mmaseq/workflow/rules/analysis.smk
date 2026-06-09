# Assemblers

rule assembly:
    input:
        input_assembly = "%s/{sample}/raw/{assembler}/{sample}.fasta" %outdir
    output:
        output_assembly = "%s/{sample}/Assemblies/{sample}_{assembler}.fasta" %outdir
    log:
        stdout = "%s/Assemblies/{sample}_{assembler}_assembly.log" %logdir
    message:
        "[assembly]: Moving {wildcards.assembler} assembly for {wildcards.sample} from raw location to assembly folder"
 
    shell:
        """
        mkdir -p $(dirname {output.output_assembly})

        cmd="cp {input.input_assembly} {output.output_assembly}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 
        """

# Gene mappers

rule assembly_minimap2:
    input:
        assembly = rules.assembly.output.output_assembly,
        database = rules.fetch_genbank.output.fasta
    params:
        options = lambda wc: sample_configs[wc.sample]["assembly_minimap2"]["options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/minimap2/{{assembler}}_{{database}}.sam")
    conda:
        ENVS_DIR / "minimap2.yaml"
    log:
        stdout = "%s/{sample}/minimap2/{assembler}_{database}.log" %logdir
    message:
        "[assembly_minimap2]: Running Minimap2 for {wildcards.database} on {wildcards.assembler} for {wildcards.sample}"
    shell:
        r"""
        mkdir -p $(dirname {output.results})

        cmd="minimap2 {params.options} {input.database} {input.assembly} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule blastn:
    input:
        # A complete access to the wildcard is needed, if we try to call the output of different rule we have the blending of wildcards 
        assembly = rules.assembly.output.output_assembly,
        database = rules.fetch_blast_database.output.source
    params:
        options = lambda wc: sample_configs[wc.sample]["blastn"]["options"]
    output:
        results = "%s/{sample}/raw/blastn/blast_{assembler}_{database}.tsv" %outdir
    conda:
        ENVS_DIR / "blast.yaml"
    log:
        stdout = "%s/{sample}/blastn_{assembler}_{database}.log" %logdir
    message:
        "[setup_{wildcards.database}]: Setting up the {wildcards.database} database from the temporary storage folder"
    shell:
        """
        mkdir -p $(dirname {output.results})

        cmd="blastn -subject {input.database} -query {input.assembly} -outfmt '6' -out {output.results} {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule amrfinder:
    input:
        assembly = rules.assembly.output.output_assembly,
        database = rules.setup_AMRFinder.output.database
    params:
        options = lambda wc: sample_configs[wc.sample]["amrfinder"]["options"]
    output:
        result = "%s/{sample}/raw/amrfinder/{assembler}.tsv" %outdir
    conda:
        ENVS_DIR / "amrfinder.yaml"
    log:
        stdout = "%s/{sample}/amrfinder_{assembler}.log" %logdir
    message:
        "[AMRFinderPlus]: Running AMRFinderPlus for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        mkdir -p $(dirname {output.result})

        cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.options} --output {output.result}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule mlst:
    input:
        assembly = rules.assembly.output.output_assembly
    output:
        mlst_file = "%s/{sample}/raw/mlst/{assembler}_mlst.tsv" %outdir,
        mlst_tmp = temp("%s/{sample}/raw/mlst/{assembler}_mlst.mp" %outdir)
    conda:
        ENVS_DIR / "mlst.yaml"
    log:
        stdout = "%s/{sample}/{assembler}_mlst.log" %logdir
    message:
        "[MLST]: Running MLST on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.mlst_file})

        cmd="mlst {input.assembly} --label $(basename {input.assembly} .fasta) > {output.mlst_tmp}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        awk -f {SCRIPTS_DIR}/mlst_header.awk {output.mlst_tmp} > {output.mlst_file}
        """


rule LREfinder:
    input:
        res = rules.kmeraligner.output.results,
        matrix = rules.kmeraligner.output.matrix
    # params:
    #     options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["blastn"]["options"],    
    output:
        results = "%s/{sample}/raw/LREfinder/{database}.tsv" %outdir,
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/LRE-finder_{database}.log" %logdir
    message:
        "[LRE-finder]: Identify genes and mutations leading to linezolid resistance in E. faecalis and E. faecium"
    shell:
        """
        mkdir -p $(dirname {output.results})
    
        cmd="python {SCRIPTS_DIR}/LRE-Typer.py -ires {input.res} -imat {input.matrix} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "LRE-finder successfully executed" > {log.stdout}
        """


rule chtyper:
    input:
        results = rules.kmeraligner.output.results
    params:
        id = 90,
        coverage = 60
    output:
        filtered_tsv = "%s/{sample}/raw/chtyper/{database}_chtyper.tsv" % outdir
    log:
        stdout = "%s/{sample}/{database}_chtyper.log" %logdir
    message:
        "[CH Typer]: Running Chtyper on {wildcards.database} assembly for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        echo "Running awk filter on {input.results}" > {log.stdout} 2>&1

        awk -F'\t' 'NR==1{{for(i=1;i<=NF;i++){{if($i=="Template_Identity")id=i;if($i=="Template_Coverage")cov=i}}print;next}} ($id+0>{params.id} && $cov+0>{params.coverage})' {input.results} > {output.filtered_tsv} 2>> {log.stdout}
        """

# Species specific mappers

rule kleborate:
    input:
        assembly = rules.assembly.output.output_assembly,
        version_db = rules.setup_kleborate_amrfinder.output.version_db
    output:
        kleborate = "%s/{sample}/raw/kleborate/{assembler}/Kleborate_long.tsv" %outdir
    params:
        options = lambda wildcards: sample_configs[wildcards.sample]["kleborate"]["options"]
    conda:
        ENVS_DIR / "kleborate.yaml"
    log:
        stdout = "%s/{sample}/Kleborate_{assembler}.log" %logdir
    message:
        "[Kleborate]: Running Kleborate on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.kleborate})
        #mkdir -p $outdir

        cmd="kleborate --assemblies {input.assembly} --outdir $outdir {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # Creating long table
        (
            echo -e "Sample\tModule\tFile\tRow\tColumn\tValue"
            for f in $outdir/*.txt; do
                fname=$(basename "$f" .tsv)

                awk -v file="$fname" 'BEGIN{{FS=OFS="\t"}}
                NR==1 {{
                    for (i=2; i<=NF; i++) header[i]=$i
                    next
                }}
                {{
                    rownum++
                    for (i=2; i<=NF; i++) {{
                        print $1, "kleborate", file, rownum, header[i], $i
                    }}
                }}' "$f"
            done
        ) > {output.kleborate}
        """


rule meningotype:
    input:
        assembly = rules.assembly.output.output_assembly,
    output:
        meningotype = "%s/{sample}/raw/meningotype/{assembler}_meningotype.tsv" %outdir
    conda:
        ENVS_DIR / "meningotype.yaml"
    log:
        stdout = "%s/{sample}/{assembler}_meningotype.log" %logdir
    message:
        "[Meningotype]: Running Meningotype on {wildcards.assembler} assembly for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.meningotype})

        cmd="meningotype --all {input.assembly} > {output.meningotype}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule sistr:
    input:
        assembly = rules.assembly.output.output_assembly,
        serovarlist = rules.fetch_Senterica_Serovar.output.source
    output:
        sistr_tab = "%s/{sample}/raw/sistr/{assembler}_sistr.tab" %outdir,
        gmlst_profile = "%s/{sample}/raw/sistr/{assembler}_cgmlst_profiles.csv" %outdir,
        allele_results = "%s/{sample}/raw/sistr/{assembler}_allele-results.json" %outdir
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 1
    conda:
        ENVS_DIR / "sistr.yaml"
    log:
        stdout = "%s/{sample}/{assembler}_SISTR_serovar.log" %logdir
    message:
        "[Salmonella_serovar]: Predict Salmonella serovar with SISTR"
    shell:
        """
        mkdir -p $(dirname {output.sistr_tab})

        cmd="sistr -f tab --qc -t {threads} -l {input.serovarlist} --cgmlst-profiles {output.gmlst_profile} --alleles-output {output.allele_results} --output-prediction {output.sistr_tab} {input.assembly}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule spa_typing:
    input:
        assembly = rules.assembly.output.output_assembly,
        database = rules.setup_Spatyper.output.database
    output:
        spatyper = "%s/{sample}/raw/spatyper/{assembler}_spatype_results.tsv" %outdir
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/spatyper_{assembler}.log" %logdir
    message:
        "[Spatyping]: Running Spatyper for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        outdir=$(dirname {output.spatyper})
        cmd="python {SCRIPTS_DIR}/SPATyper_V2.py -a {input.assembly} -d {input.database} -o {output.spatyper} -b $outdir/seq_db -l $outdir/spatyper.log "

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# SAM analysis

rule samtools_sam_filtration:
    input:
        sam = "%s/{sample}/raw/samtools/{database}.sam" %outdir
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["view_options"]
    output:
        bam = temp("%s/{sample}/raw/samtools/{database}_filtered.bam" %outdir)
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = "%s/{sample}/samtools_filtration_{database}.log" %logdir
    message:
        "[samtools_sam_filtration]: Filtering SAM of {wildcards.sample} mapped against {wildcards.database}"
    shell:
        """
        cmd="samtools view {input.sam} {params.options} -F 4 -bo {output.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule samtools_bam_filtration:
    input:
        bam = "%s/{sample}/raw/samtools/{database}.bam" %outdir
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["view_options"]
    output:
        bam = temp("%s/{sample}/raw/samtools/{database}_filtered.bam" %outdir)
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = "%s/{sample}/samtools_filtration_{database}.log" %logdir
    message:
        "[samtools_bam_filtration]: Filtering BAM of {wildcards.sample} mapped against {wildcards.database}"
    shell:
        """
        cmd="samtools view {input.bam} {params.options} -F 4 -bo {output.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule samtools_sort:
    input:
        bam = "%s/{sample}/raw/samtools/{database}_filtered.bam" %outdir
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["sort_options"]
    output:
        bam_sort = temp("%s/{sample}/raw/samtools/{database}_sorted.bam" %outdir),
        index = temp("%s/{sample}/raw/samtools/{database}_sorted.bam.bai" %outdir)
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = "%s/{sample}/samtools_sort_{database}.log" %logdir
    message:
        "[samtools_sort]: Sorting filtered bam for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools sort -o {output.bam_sort} {input.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="samtools index {output.bam_sort}"

        echo "\nIndexing Bam:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_pileup:
    input:
        bam_sort = rules.samtools_sort.output.bam_sort,
        reference = "%s/samtools/{database}.fasta" %database_dir
    output:
        pileup = temp("%s/{sample}/raw/bcftools/{database}_pileup.bcf" %outdir),
        index = temp("%s/{sample}/raw/bcftools/{database}_pileup.bcf.csi" %outdir)
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = "%s/{sample}/bcftools_pileup_{database}.log" %logdir
    message:
        "[bcftools_pileup]: Generating mpileup for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools mpileup -Ob -f {input.reference} {input.bam_sort} -o {output.pileup}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.pileup} -o {output.index}"

        echo "\nIndexing Pileup:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_filter_indels:
    input:
        pileup = rules.bcftools_pileup.output.pileup,
        pileup_index = rules.bcftools_pileup.output.index,
    params:
        options = lambda wc: sample_configs[wc.sample]["bcftools"]["view_options"]
    output:
        indels = temp("%s/{sample}/raw/bcftools/{database}_pileup_indels.bcf" %outdir),
        index = temp("%s/{sample}/raw/bcftools/{database}_pileup_indels.bcf.csi" %outdir)
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = "%s/{sample}/bcftools_filter_indels_{database}.log" %logdir
    message:
        "[bcftools_filter_indels]: Filtering indels of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools view {params.options} -Ob -o {output.indels} {input.pileup}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.indels} -o {output.index}"

        echo "\nIndexing Pileup:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_variant_call:
    input:
        pileup = rules.bcftools_pileup.output.pileup,
        pileup_index = rules.bcftools_pileup.output.index,
    output: 
        variants = temp("%s/{sample}/raw/bcftools/{database}_call_variants.bcf" %outdir),
        index = temp("%s/{sample}/raw/bcftools/{database}_call_variants.bcf.csi" %outdir)
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = "%s/{sample}/bcftools_variant_call_{database}.log" %logdir
    message:
        "[bcftools_variant_call]: Calling variant of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools call -mv -Ob --ploidy 1 {input.pileup} -o {output.variants}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.variants} -o {output.index}"

        echo "\nIndexing Call:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

# Custom analysis
rule cdiff_repeat_identifier:
    input:
        seqs  = expand(rules.fetch_type_repeat_sequence.output.seq, TR = ["TR6", "TR10"]),
        metas = expand(rules.fetch_type_repeat_metadata.output.meta, TR = ["TR6", "TR10", "TRST"]),
        assembly = rules.assembly.output.output_assembly
    params:
        repeats = lambda wc: sample_configs[wc.sample]["cdiff_repeat_identifier"]["repeats"],
        combos = lambda wc: sample_configs[wc.sample]["cdiff_repeat_identifier"]["combos"]
    output:
        repeat_types = "%s/{sample}/raw/cdiff_repeat_identifier/{assembler}_repeat_types.tsv" %outdir
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/cdiff_repeat_identifier_{assembler}_repeat_types.log" %logdir
    message:
        "[CDiff Repeat identifier]: Identifying C. Difficile repeats in {wildcards.sample} on {wildcards.assembler} assembly"
    shell:
        """
        mkdir -p $(dirname {output.repeat_types})

        db_dir=$(dirname {input.seqs} | uniq)

        cmd="python {SCRIPTS_DIR}/Repeat_Identifier.py --fasta {input.assembly} --ref_seq {input.seqs} --ref_meta {input.metas} --output {output.repeat_types} --sample_id {wildcards.sample} --repeats {params.repeats} --combos {params.combos} --suffix tsv"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 
        """


rule snp_identifier:
    input:
        variants = rules.bcftools_variant_call.output.variants,
        variants_index = rules.bcftools_variant_call.output.index,
    params:
        options = lambda wc: sample_configs[wc.sample]["snp_identifier"]["options"],
        metafile = "%s/SNP_metafile.tsv" %SCREENING_DIR
    output:
        indentified_variants = "%s/{sample}/raw/snp_identifier/{database}.tsv" %outdir
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/snp_identifier_{database}.log" %logdir
    message:
        "[SNP Identifier]: Identifying SNPs of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="python {SCRIPTS_DIR}/SNP_identifier.py {params.options} --call {input.variants} --metafile {params.metafile} --output {output.indentified_variants}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule deletion_identifier:
    input:
        kma_seq = rules.kmeraligner_consensus.output.seq,
        indels = rules.bcftools_filter_indels.output.indels,
        indels_index = rules.bcftools_filter_indels.output.index,
        variants = rules.bcftools_variant_call.output.variants,
        variants_index = rules.bcftools_variant_call.output.index,
        asm_aln = rules.assembly_minimap2.output.results
    params:
        options  = lambda wc: sample_configs[wc.sample]["deletion_identifier"]["options"],
        metafile = "%s/deletion_metafile.tsv" %SCREENING_DIR
    output:
        identified_variants = f"{outdir}/{{sample}}/raw/deletion_identifier/{{assembler,[^_]+}}_{{database}}.tsv" #added regex expression to ensure assemblies cannot contain '_' which our database also does
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/deletion_identifier_{assembler}_{database}.log" %logdir
    message:
        "[Deletion Identifier]: Identifying deletions of {wildcards.database} on {wildcards.sample} ({wildcards.assembler})"
    shell:
        """
        cmd="python {SCRIPTS_DIR}/deletion_identifier.py {params.options} --fsa {input.kma_seq} --call {input.variants} --mpileup {input.indels} --metafile {params.metafile} --sam {input.asm_aln} --output {output.identified_variants}"


        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


# Custom data wranglers

rule kma_filter:
    input:
        results = rules.kmeraligner.output.results,
        database = rules.setup_kmeraligner_index.output.names
    params:
        options = lambda wildcards: sample_configs[wildcards.sample]["kma_filter"]["options"],
        metafile = "%s/kma_filter.tsv" %SCREENING_DIR
    output:
        filtered_tsv = "%s/{sample}/raw/kma_filter/{database}.tsv" % outdir
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = "%s/{sample}/KMA_results/{sample}_{database}.log" %logdir
    message:
        "[kma_filter]: Filtering KMA .res result for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        cmd="python {SCRIPTS_DIR}/KMA_Filter.py --KMA_res {input.results} --metafile {params.metafile} --sample_id {wildcards.sample} --output {output.filtered_tsv} {params.options} > {log.stdout} 2>&1"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """