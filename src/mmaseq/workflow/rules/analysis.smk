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


### Read mapping tools ###
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
        results = f"{outdir}/{{sample}}/raw/amrfinder/amrfinder_{{assembler}}.tsv"
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
        
        cmd="amrfinder --nucleotide {input.assembly} --database {input.database} {params.options} --output {output.results}"

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
    params:
        tmp_results = f"sistr_{{assembler}}.tsv.tab"
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


rule seqsero2:
    input:
        assembly = rules.assembly.output.assembly
    params:
        tmp_results = "SeqSero_result.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/seqsero2/seqsero2_{{assembler}}.tsv"
    threads: max(1, workflow.cores - 1 - (workflow.cores - 1) % 2)
    priority: 1
    conda:
        ENVS_DIR / "seqsero2.yaml"
    log:
        stdout = f"{logdir}/seqsero2_{{assembler}}_{{sample}}.log"
    message:
        "[seqsero2]: Running seqsero2 on {wildcards.sample} using {wildcards.assembler} contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        mkdir -p $OUTDIR

        cmd="SeqSero2_package.py -m k -t 4 -b mem -i {input.assembly} -d $OUTDIR -n {wildcards.sample} -p {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule plasmidfinder:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.setup_plasmidfinder.output.database
    params:
        tmp_results = "results_tab.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/plasmidfinder/plasmidfinder_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "plasmidfinder.yaml"
    log:
        stdout = f"{logdir}/plasmidfinder_{{assembler}}_{{sample}}.log"
    message:
        "[plasmidfinder]: Running PlasmidFinder on {wildcards.sample} using {wildcards.assembler} contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})

        cmd="plasmidfinder.py -i {input.assembly} -o $OUTDIR -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1     

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule resfinder:
    input:
        assembly = rules.assembly.output.assembly,
        res_database = rules.setup_resfinder.output.database
    params:
        tmp_results = "ResFinder_results_tab.txt",
        options = lambda wc: sample_configs[wc.sample]["resfinder"]["options"]
    output:
        results = f"{outdir}/{{sample}}/raw/SR/resfinder/resfinder_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = f"{logdir}/resfinder_{{assembler}}_{{sample}}.log"
    message:
        "[resfinder]: Running ResFinder on {wildcards.sample} using {wildcards.assembler} contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        cmd="run_resfinder.py -ifa {input.assembly} -o $OUTDIR --acquired -db_res {input.res_database} {params.options}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule pointfinder:
    input:
        assembly = rules.assembly.output.assembly,
        res_database = rules.setup_resfinder.output.database,
        point_database = rules.setup_pointfinder.output.database
    params:
        tmp_results = "PointFinder_results.txt",
        options = lambda wc: sample_configs[wc.sample]["pointfinder"]["options"]
    output:
        results = f"{outdir}/{{sample}}/raw/SR/pointfinder/pointfinder_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = f"{logdir}/pointfinder_{{assembler}}_{{sample}}.log"
    message:
        "[pointfinder]: Running PointFinder on {wildcards.sample} using {wildcards.assembler} contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        cmd="run_resfinder.py -ifa {input.assembly} -o $OUTDIR -db_res {input.res_database} --point -db_point {input.point_database} {params.options}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule disinfinder:
    input:
        assembly = rules.assembly.output.assembly,
        res_database = rules.setup_resfinder.output.database,
        disin_database = rules.setup_disinfinder.output.database
    params:
        tmp_results = "DisinFinder_results_tab.txt",
        options = lambda wc: sample_configs[wc.sample]["disinfinder"]["options"]
    output:
        results = f"{outdir}/{{sample}}/raw/SR/disinfinder/disinfinder_{{assembler}}.tsv"
    conda:
        ENVS_DIR / "resfinder.yaml"
    log:
        stdout = f"{logdir}/disinfinder_{{assembler}}_{{sample}}.log"
    message:
        "[disinfinder]: Running DisinFinder on {wildcards.sample} using {wildcards.assembler} contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        cmd="run_resfinder.py -ifa {input.assembly} -o $OUTDIR -db_res {input.res_database} --disinfectant -db_disinf {input.disin_database} {params.options}"
    
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule virulencefinder:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.setup_virulencefinder.output.database
    params:
        tmp_results = "results_tab.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/SR/virulencefinder/virulencefinder_{{assembler}}.tsv",
    conda:
        ENVS_DIR / "virulencefinder.yaml"
    log:
        stdout = f"{logdir}/virulencefinder_{{assembler}}_{{sample}}.log"
    message:
        "[virulencefinder]: Running VirulenceFinder on {wildcards.sample} using {wildcards.assembler} contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})
        cmd="python -m virulencefinder -ifa {input.assembly} -o $OUTDIR -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """


rule serotypefinder:
    input:
        assembly = rules.assembly.output.assembly,
        database = rules.setup_serotypefinder.output.database
    params:
        tmp_results = "results_tab.tsv"
    output:
        results = f"{outdir}/{{sample}}/raw/SR/serotypefinder/serotypefinder_{{assembler}}.tsv",
    conda:
        ENVS_DIR / "serotypefinder.yaml"
    log:
        stdout = f"{logdir}/serotypefinder_{{assembler}}_{{sample}}.log"
    message:
        "[serotypefinder]: Running SerotypeFinder on {wildcards.sample} using {wildcards.assembler} contigs"
    shell:
        """
        OUTDIR=$(dirname {output.results})

        cmd="serotypefinder -i {input.R1} -o $OUTDIR -p {input.database} -x"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Renaming result files" >> {log.stdout} 2>&1
        mv $OUTDIR/{params.tmp_results} {output.results} >> {log.stdout} 2>&1
        """