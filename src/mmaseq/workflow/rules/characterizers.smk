# Rule MLST:
# Runs Multi Locus Sequence Type to determine the ST profile of isolate
rule mlst:
    input:
        assembly = rules.assembly.output.output_assembly
    output:
        mlst_file = "%s/{sample}/mlst/{assembler}_mlst.tsv" %outdir,
        mlst_tmp = temp("%s/{sample}/mlst/{assembler}_mlst.mp" %outdir),
        tool_version = "%s/{sample}/mlst/{assembler}_mlst_version.txt" %outdir,
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

        # 2) create version file with date
        version_cmd="mlst --version"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
    	"""

# Rule Kleborate:
# Runs Kleborate characterising virulence and resistance in pathogen assemblies
rule kleborate:
    input:
        assembly = rules.assembly.output.output_assembly,
        version_db = rules.setup_kleborate_amrfinder.output.version_db
    output:
        outfile = helper_functions.determine_rule_output(outdir, "{sample}", "kleborate", results_catalogue, sample_configs)
    params:
        options = lambda wildcards: sample_configs[wildcards.sample]["kleborate"]["options"],
        output = lambda wildcards: sample_configs[wildcards.sample]["kleborate"]["output"]
    conda:
        ENVS_DIR / "kleborate.yaml"
    log:
    	stdout = "%s/{sample}/Kleborate_{assembler}.log" %logdir
    message:
    	"[Kleborate]: Running Kleborate on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        #mkdir -p $outdir

        cmd="kleborate --assemblies {input.assembly} --outdir $(dirname {output.kleborate}) {params.options} && mv {params.output} {output.outfile}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    	"""

rule chtyper:
    input:
        results = rules.custom_kmeralignment.output.results
    params:
        id = 90,
        coverage = 60
    output:
        filtered_tsv = "%s/{sample}/chtyper/{database}_chtyper.tsv" % outdir
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


# Rule meningotype
# Runs meningotype on SPAdes assembled contigs to perform serotyping of N. Meningmeningitidis contigs
rule meningotype:
    input:
        assembly = rules.assembly.output.output_assembly,
    output:
        meningotype = "%s/{sample}/meningotype/{assembler}_meningotype.tsv" %outdir,
        tool_version = "%s/{sample}/meningotype/{assembler}_meningotype_version.txt" %outdir,
    conda:
        ENVS_DIR / "meningotype.yaml"
    log:
        stdout = "%s/{sample}/{assembler}_meningotype.log" %logdir
    message:
    	"[Meningotype]: Running Meningotype on {wildcards.assembler} assembly for {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.meningotype})

        cmd="meningotype --all {input.assembly} > {output.meningotype} 2> {log.stdout}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # 2) create version file with date
        version_cmd="meningotype --version"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
    	"""


rule seqsero2:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        seqsero = "%s/{sample}/seqsero2/SeqSero_result.tsv" %outdir,
        tool_version = "%s/{sample}/seqsero2/SeqSero_version.txt" %outdir,
    threads: workflow.cores - 1 - (workflow.cores - 1) % 2
    priority: 1
    conda:
        ENVS_DIR / "seqsero2.yaml"
    log:
        stdout = "%s/{sample}/seqsero2.log" %logdir
    message:
        "[seqsero2]: Running seqsero2 on {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.seqsero})
        mkdir -p $outdir

        cmd="SeqSero2_package.py -m k -t 2 -b mem -i {input.R1} {input.R2} -d $outdir -n {wildcards.sample} -p {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # 2) create version file with date
        version_cmd="SeqSero2_package.py -v"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
        """

rule sistr:
    input:
        assembly = rules.assembly.output.output_assembly,
        serovarlist = rules.fetch_Senterica_Serovar.output.source
    output:
        sistr_tab = "%s/{sample}/sistr/{assembler}_sistr.tab" %outdir,
        gmlst_profile = "%s/{sample}/sistr/{assembler}_cgmlst_profiles.csv" %outdir,
        allele_results = "%s/{sample}/sistr/{assembler}_allele-results.json" %outdir,
        tool_version = "%s/{sample}/sistr/{assembler}_sistr_version.txt" %outdir,
    threads: workflow.cores - 1 - (workflow.cores - 1) % 2
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
        
        # 2) create version file with date
        version_cmd="sistr --version 2>/dev/null"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
        """
