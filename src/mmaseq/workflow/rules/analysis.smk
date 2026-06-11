### Assembly ###

rule assembly:
    input:
        assembly = f"{outdir}/{{sample}}/raw/{{assembler}}/{{assembler}}_{{sample}}.fasta"
    output:
        assembly = f"{outdir}/{{sample}}/Assemblies/{{assembler}}_{{sample}}.fasta"
    log:
        stdout = f"{logdir}/sync_{{assembler}}_{{sample}}.log"
    message:
        "[assembly]: Syncronizing {wildcards.assembler} for {wildcards.sample} from raw location to assembly folder"
 
    shell:
        """
        mkdir -p $(dirname {output.assembly})

        cmd="rsync -av {input.assembly} {output.assembly}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 
        """

### Mapping ###

rule blastn:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.fetch_custom_blast_database.output.source
    params:
        options = lambda wc: sample_configs[wc.sample]["blastn"]["options"]
    output:
        results = f"{outdir}/{{sample}}/raw/blastn/blastn_{{database}}_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "blast.yaml"
    log:
        stdout = f"{logdir}/blastn_{{database}}_{{assembler}}_{{sample}}.log"
    message:
        "[blastn]: Blasting {wildcards.database} against {wildcards.sample} database from the temporary storage folder"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="blastn -subject {input.database} -query {input.assembly} -outfmt '6' -out {output.results} {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule minimap2:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.fetch_genbank.output.fasta
    params:
        options = lambda wc: sample_configs[wc.sample]["minimap2"]["options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/minimap2/minimap2_{{assembler}}_{{database}}.sam")
    conda:
        ENVS_DIR / "minimap2.yaml"
    log:
        stdout = f"{logdir}/minimap2_{{assembler}}_{{database}}_{{sample}}.log"
    message:
        "[minimap2]: Running Minimap2 for {wildcards.database} on {wildcards.assembler} for {wildcards.sample}"
    shell:
        r"""
        mkdir -p $(dirname {output.results})

        cmd="minimap2 {params.options} {input.database} {input.assembly} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

### Read mapping tools ###

rule samtools_sam_filtration:
    input:
        sam = f"{outdir}/{{sample}}/raw/samtools/{{database}}.sam"
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["view_options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/samtools/samtools_bam_filtration_{{database}}.bam")
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = f"{logdir}/samtools_sam_filtration_{{database}}_{{sample}}.log"
    message:
        "[samtools_sam_filtration]: Filtering kmeralignment output for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools view {input.sam} {params.options} -F 4 -bo {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule samtools_bam_filtration:
    input:
        bam = f"{outdir}/{{sample}}/raw/samtools/{{database}}.bam"
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["view_options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/samtools/samtools_bam_filtration_{{database}}.bam")
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = f"{logdir}/samtools_bam_filtration_{{database}}_{{sample}}.log"
    message:
        "[samtools_bam_filtration]: Filtering kmeralignment output for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools view {input.bam} {params.options} -F 4 -bo {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule samtools_sort:
    input:
        bam = f"{outdir}/{{sample}}/raw/samtools/samtools_bam_filtration_{{database}}.bam"
    params:
        options = lambda wc: sample_configs[wc.sample]["samtools"]["sort_options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/samtools/samtools_sort_{{database}}.bam"),
        index = temp(f"{outdir}/{{sample}}/raw/samtools/samtools_sort_{{database}}.bam.bai")
    conda:
        ENVS_DIR / "samtools.yaml"
    log:
        stdout = f"{logdir}/samtools_sort_{{database}}_{{sample}}.log"
    message:
        "[samtools_sort]: Sorting filtered bam for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="samtools sort -o {output.results} {input.bam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="samtools index {output.results}"

        echo "\nIndexing Bam:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_pileup:
    input:
        bam_sort = rules.samtools_sort.output.results,
        reference = f"{database_dir}/samtools/{{database}}.fasta"
    output:
        results = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_pileup_{{database}}.bcf"),
        index = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_pileup_{{database}}.bcf.csi")
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = f"{logdir}/bcftools_pileup_{{database}}_{{sample}}.log"
    message:
        "[bcftools_pileup]: Generating mpileup for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools mpileup -Ob -f {input.reference} {input.bam_sort} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.results} -o {output.index}"

        echo "\nIndexing Pileup:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_filter_indels:
    input:
        pileup = rules.bcftools_pileup.output.results,
        pileup_index = rules.bcftools_pileup.output.index,
    params:
        options = lambda wc: sample_configs[wc.sample]["bcftools"]["view_options"]
    output:
        results = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_filter_indels_{{database}}.bcf"),
        index = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_filter_indels_{{database}}.bcf.csi")
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = f"{logdir}/bcftools_filter_indels_{{database}}_{{sample}}.log"
    message:
        "[bcftools_filter_indels]: Filtering indels of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools view {params.options} -Ob -o {output.results} {input.pileup}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.results} -o {output.index}"

        echo "\nIndexing Pileup:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule bcftools_variant_call:
    input:
        pileup = rules.bcftools_pileup.output.results,
        pileup_index = rules.bcftools_pileup.output.index,
    output: 
        results = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_variant_call_{{database}}.bcf"),
        index = temp(f"{outdir}/{{sample}}/raw/bcftools/bcftools_variant_call_{{database}}.bcf.csi")
    conda:
        ENVS_DIR / "bcftools.yaml"
    log:
        stdout = f"{logdir}/bcftools_variant_call_{{database}}_{{sample}}.log"
    message:
        "[bcftools_variant_call]: Calling variant of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="bcftools call -mv -Ob --ploidy 1 {input.pileup} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="bcftools index -f {output.results} -o {output.index}"

        echo "\nIndexing Call:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """

### SNP analysis ###

rule snp_identifier:
    input:
        variants = rules.bcftools_variant_call.output.results,
        variants_index = rules.bcftools_variant_call.output.index,
    params:
        options = lambda wc: sample_configs[wc.sample]["snp_identifier"]["options"],
        metafile = f"{SCREENING_DIR}/SNP_metafile.tsv"
    output:
        indentified_variants = f"{outdir}/{{sample}}/raw/snp_identifier/snp_identifier_{{database}}.tsv"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/snp_identifier_{{database}}_{{sample}}.log"
    message:
        "[snp_identifier]: Identifying SNPs of {wildcards.database} on {wildcards.sample}"
    shell:
        """
        cmd="python {SCRIPTS_DIR}/SNP_identifier.py {params.options} --call {input.variants} --metafile {params.metafile} --output {output.indentified_variants}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule deletion_identifier:
    input:
        consensus_seq = rules.kmeraligner_consensus.output.results,
        indels = rules.bcftools_filter_indels.output.results,
        indels_index = rules.bcftools_filter_indels.output.index,
        variants = rules.bcftools_variant_call.output.results,
        variants_index = rules.bcftools_variant_call.output.index,
        asm_aln = rules.minimap2.output.results
    params:
        options  = lambda wc: sample_configs[wc.sample]["deletion_identifier"]["options"],
        metafile = f"{SCREENING_DIR}/deletion_metafile.tsv"
    output:
        identified_variants = f"{outdir}/{{sample}}/raw/deletion_identifier/deletion_identifier_{{database}}_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/deletion_identifier_{{database}}_{{assembler}}_{{sample}}.log"
    message:
        "[deletion_identifier]: Identifying deletions of {wildcards.database} on {wildcards.sample} ({wildcards.assembler})"
    shell:
        """
        cmd="python {SCRIPTS_DIR}/deletion_identifier.py {params.options} --fsa {input.consensus_seq} --call {input.variants} --mpileup {input.indels} --metafile {params.metafile} --sam {input.asm_aln} --output {output.identified_variants}"


        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule cdiff_repeat_identifier:
    input:
        seqs  = expand(rules.fetch_type_repeat_sequence.output.seq, TR = ["TR6", "TR10"]),
        metas = expand(rules.fetch_type_repeat_metadata.output.meta, TR = ["TR6", "TR10", "TRST"]),
        assembly = rules.assembly.output.assembly
    params:
        repeats = lambda wc: sample_configs[wc.sample]["cdiff_repeat_identifier"]["repeats"],
        combos = lambda wc: sample_configs[wc.sample]["cdiff_repeat_identifier"]["combos"]
    output:
        repeat_types = f"{outdir}/{{sample}}/raw/cdiff_repeat_identifier/cdiff_repeat_identifier_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/cdiff_repeat_identifier_{{assembler}}_{{sample}}.log"
    message:
        "[cdiff_repeat_identifier]: Identifying C. Difficile repeats in {wildcards.sample} on {wildcards.assembler} assembly"
    shell:
        """
        mkdir -p $(dirname {output.repeat_types})

        db_dir=$(dirname {input.seqs} | uniq)

        cmd="python {SCRIPTS_DIR}/Repeat_Identifier.py --fasta {input.assembly} --ref_seq {input.seqs} --ref_meta {input.metas} --output {output.repeat_types} --sample_id {wildcards.sample} --repeats {params.repeats} --combos {params.combos} --suffix tsv"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 
        """


### Characterizers ###

rule amrfinder:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.setup_amrfinder.output.database
    params:
        options = lambda wc: sample_configs[wc.sample]["amrfinder"]["options"]
    output:
        result = f"{outdir}/{{sample}}/raw/amrfinder/amrfinder_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "amrfinder.yaml"
    log:
        stdout = f"{logdir}/amrfinder_{{assembler}}_{{sample}}.log"
    message:
        "[amrfinder]: Running AMRFinderPlus for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR
        
        cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.options} --output {output.result}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule lrefinder:
    input:
        res = rules.kmeraligner.output.results,
        matrix = rules.kmeraligner.output.matrix
    output:
        results = f"{outdir}/{{sample}}/raw/lrefinder/lrefinder_{{database}}.tsv",
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/lrefinder_{{database}}_{{sample}}.log"
    message:
        "[LREfinder]: Identify genes and mutations for linezolid resistance in {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="python {SCRIPTS_DIR}/LRE-Typer.py -ires {input.res} -imat {input.matrix} -o {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule mlst:
    input:
        assembly = rules.assembly.output.assembly
    output:
        results = f"{outdir}/{{sample}}/raw/mlst/mlst_{{assembler}}.tsv",
        mlst_tmp = temp(f"%s/{{sample}}/raw/mlst/mlst_{{assembler}}.tmp")
    conda:
        ENVS_DIR / "mlst.yaml"
    log:
    	stdout = f"{logdir}/mlst_{{assembler}}_{{sample}}.log"
    message:
    	"[mlst]: Running MLST on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.results})

        cmd="mlst {input.assembly} --label $(basename {input.assembly} .fasta) > {output.mlst_tmp}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        awk -f {SCRIPTS_DIR}/mlst_header.awk {output.mlst_tmp} > {output.results}
    	"""


rule kleborate:
    input:
        assembly = rules.assembly.output.assembly,
        version_db = rules.setup_kleborate_amrfinder.output.version_db
    output:
        results = f"{outdir}/{{sample}}/raw/kleborate/kleborate_{{assembler}}.tsv"
    params:
        options = lambda wildcards: sample_configs[wildcards.sample]["kleborate"]["options"]
    conda:
        ENVS_DIR / "kleborate.yaml"
    log:
    	stdout = f"{logdir}/kleborate_{{assembler}}_{{sample}}.log"
    message:
    	"[kleborate]: Running Kleborate on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        #mkdir -p $OUTDIR

        cmd="kleborate --assemblies {input.assembly} --outdir $OUTDIR {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # Creating long table
        (
            echo -e "Sample\tModule\tFile\tRow\tColumn\tValue"
            for f in $OUTDIR/*.txt; do
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
        ) > {output.results}
    	"""


rule chtyper:
    input:
        results = rules.kmeraligner.output.results
    params:
        id = 90,
        coverage = 60
    output:
        results = f"{outdir}/{{sample}}/raw/chtyper/chtyper_{{database}}.tsv"
    log:
        stdout = f"{logdir}/chtyper_{{database}}_{{sample}}.log"
    message:
    	"[chtyper]: Running Chtyper on {wildcards.database} assembly for {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        echo "Running awk filter on {input.results}" > {log.stdout} 2>&1

        awk -F'\t' 'NR==1{{for(i=1;i<=NF;i++){{if($i=="Template_Identity")id=i;if($i=="Template_Coverage")cov=i}}print;next}} ($id+0>{params.id} && $cov+0>{params.coverage})' {input.results} > {output.results} 2>> {log.stdout}
        """


rule meningotype:
    input:
        assembly = rules.assembly.output.assembly,
    output:
        results = f"{outdir}/{{sample}}/raw/meningotype/meningotype_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "meningotype.yaml"
    log:
        stdout = f"{logdir}/meningotype_{{assembler}}_{{sample}}.log"
    message:
    	"[meningotype]: Running Meningotype on {wildcards.assembler} assembly for {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="meningotype --all {input.assembly} > {output.results}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    	"""


rule sistr:
    input:
        assembly = rules.assembly.output.assembly,
        serovarlist = rules.fetch_senterica_serovar.output.source
    output:
        results = f"{outdir}/{{sample}}/raw/sistr/sistr_{{assembler}}.tsv",
        cgmlst = temp(f"{outdir}/{{sample}}/raw/sistr/sistr_cgmlst_profiles_{{assembler}}.csv"), # Not listed anywhere else, kept for history sake...
        alleles = temp(f"{outdir}/{{sample}}/raw/sistr/sistr_alleles_{{assembler}}.json") # Not listed anywhere else...
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 1
    conda:
        ENVS_DIR / "sistr.yaml"
    log:
        stdout = f"{logdir}/sistr_{{assembler}}_{{sample}}.log"
    message:
        "[sistr]: Predict Salmonella serovar with SISTR"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="sistr -f tab --qc -t {threads} -l {input.serovarlist} --cgmlst-profiles {output.cgmlst} --alleles-output {output.alleles} --output-prediction {output.results} {input.assembly}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


rule spa_typing:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.setup_spatyper.output.database
    output:
        results = f"{outdir}/{{sample}}/raw/spatyper/spa_typing_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/spa_typing_{{assembler}}_{{sample}}.log"
    message:
        "[spa_typing]: Running Spatyper for {wildcards.sample} using ({wildcards.assembler}) contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        cmd="python {SCRIPTS_DIR}/SPATyper_V2.py -a {input.assembly} -d {input.database} -o {output.results} -b $OUTDIR/seq_db -l $OUTDIR/spatyper.log "

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
        """


### Wranglers ###

rule kmeraligner_wrangler:
    input:
        results = rules.kmeraligner.output.results,
        database = rules.setup_kmeraligner_index.output.names
    params:
        options = lambda wildcards: sample_configs[wildcards.sample]["kmeraligner_wrangler"]["options"],
        metafile = f"{SCREENING_DIR}/kmeraligner_wrangler.tsv"
    output:
        filtered_tsv = f"{outdir}/{{sample}}/raw/kmeraligner_wrangler/kmeraligner_wrangler_{{database}}.tsv"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/kmeraligner_wrangler_{{database}}_{{sample}}.log"
    message:
        "[KMA kmeraligner_wrangler]: Filtering KMA .res result for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.filtered_tsv})

        cmd="python {SCRIPTS_DIR}/KMA_Filter.py --KMA_res {input.results} --metafile {params.metafile} --sample_id {wildcards.sample} --output {output.filtered_tsv} {params.options} > {log.stdout} 2>&1"

        echo "Executing command:\n$cmd\n" > {log.stdout}
        eval $cmd >> {log.stdout} 2>&1
        """
