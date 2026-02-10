# Rule MLST:
# Runs Multi Locus Sequence Type to determine the ST profile of isolate
rule mlst:
    input:
        assembly = rules.assembly.output.output_assembly
    output:
        mlst_file = "%s/{sample}/mlst/{assembler}_mlst.tsv" %output_folder,
        tool_version = "%s/{sample}/mlst/{assembler}_mlst_version.txt" %output_folder,
    conda:
        "../envs/mlst.yaml"
    log:
    	stdout = "Logs/{sample}/{assembler}_mlst.log"
    message:
    	"[MLST]: Running MLST on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.mlst_file})

        cmd="mlst {input.assembly} --label $(basename {input.assembly} .fasta) > {output.mlst_file}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

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
    output:
        kleborate = "%s/{sample}/kleborate/{assembler}/klebsiella_pneumo_complex_output.txt" %output_folder,
        kleborate_hAMRonization = "%s/{sample}/kleborate/{assembler}/klebsiella_pneumo_complex_hAMRonization_output.txt" %output_folder,
        tool_version = "%s/{sample}/kleborate/{assembler}/klebsiella_pneumo_version.txt" %output_folder,
    params:
        options = lambda wildcards: sample_configs[wildcards.sample]["kleborate"]["options"]
    conda:
        "../envs/kleborate.yaml"
    log:
    	stdout = "Logs/{sample}/Kleborate_{assembler}.log"
    message:
    	"[Kleborate]: Running Kleborate on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        outdir=$(dirname {output.kleborate})

        mkdir -p $outdir

        cmd="kleborate --assemblies {input.assembly} --outdir $outdir {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1
    	
        # 2) create version file with date
        version_cmd="kleborate --version 2>/dev/null"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s%s\t%s\n' "Kleborate_" "$version_str" "$date_str" > {output.tool_version}
        """

rule chtyper:
    input:
        results = rules.custom_kmeralignment.output.results
    params:
        id = 90,
        coverage = 60
    output:
        filtered_tsv = "%s/{sample}/chtyper/{database}_chtyper.tsv" % output_folder
    log:
        stdout = "Logs/{sample}/{database}_chtyper.log"
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
        meningotype = "%s/{sample}/meningotype/{assembler}_meningotype.tsv" %output_folder,
        tool_version = "%s/{sample}/meningotype/{assembler}_meningotype_version.txt" %output_folder,
    conda:
        "../envs/meningotype.yaml"
    log:
        stdout = "Logs/{sample}/{assembler}_meningotype.log"
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
        seqsero = "%s/{sample}/seqsero2/SeqSero_result.tsv" %output_folder,
        tool_version = "%s/{sample}/seqsero2/SeqSero_version.txt" %output_folder,
    threads:
        max(1, workflow.cores * 0.3333333)
    priority: 1
    conda:
        "../envs/seqsero2.yaml"
    log:
        stdout = 'Logs/{sample}/seqsero2.log'
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
        sistr_tab = "%s/{sample}/sistr/{assembler}_sistr.tab" %output_folder,
        gmlst_profile = "%s/{sample}/sistr/{assembler}_cgmlst_profiles.csv" %output_folder,
        allele_results = "%s/{sample}/sistr/{assembler}_allele-results.json" %output_folder,
        tool_version = "%s/{sample}/sistr/{assembler}_sistr_version.txt" %output_folder,
    threads:
        max(1, workflow.cores * 0.3333333)
    priority: 2
    conda:
        "../envs/sistr.yaml"
    log:
        stdout = 'Logs/{sample}/{assembler}_SISTR_serovar.log'
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
