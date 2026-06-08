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


rule custom_blaster:
    input:
        # A complete access to the wildcard is needed, if we try to call the output of different rule we have the blending of wildcards 
        assembly = rules.assembly.output.output_assembly,
        database = rules.fetch_custom_blast_database.output.source
    params:
        options = lambda wc: sample_configs[wc.sample]["custom_blaster"]["options"]
    output:
        results = "%s/{sample}/raw/custom_blaster/blast_{assembler}_{database}.tsv" %outdir
    conda:
        ENVS_DIR / "blast.yaml"
    log:
        stdout = "%s/{sample}/custom_blaster_{assembler}_{database}.log" %logdir
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
