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
        mkdir -p $(dirname {output.mlst_file})

        cmd="mlst {input.assembly} --label $(basename {input.assembly} .fasta) > {output.mlst_tmp}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        awk -f {SCRIPTS_DIR}/mlst_header.awk {output.mlst_tmp} > {output.mlst_file}
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
        results = rules.custom_kmeralignment.output.results
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


rule seqsero2:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    params:
        tmp_results = "SeqSero_result.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/seqsero2/seqsero2.tsv"
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 1
    conda:
        ENVS_DIR / "seqsero2.yaml"
    log:
        stdout = f"{logdir}/seqsero2_{{sample}}.log"
    message:
        "[seqsero2]: Running seqsero2 on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="SeqSero2_package.py -m k -t 2 -b mem -i {input.R1} {input.R2} -d $OUTDIR -n {wildcards.sample} -p {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
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


rule serotypefinder:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_serotypefinder.output.database
    params:
        tmp_results = "results_tab.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/serotypefinder/serotypefinder.tsv",
    conda:
        ENVS_DIR / "serotypefinder.yaml"
    log:
        stdout = f"{logdir}/serotypefinder_{{sample}}.log"
    message:
        "[serotypefinder]: Running SerotypeFinder on {wildcards.sample}"
    shell:
        """
        OUTDIR=$(dirname {output.results})

        cmd="serotypefinder -i {input.R1} {input.R2} -o $OUTDIR -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
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
