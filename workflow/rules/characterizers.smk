# Rule MLST:
# Runs Multi Locus Sequence Type to determine the ST profile of isolate
rule mlst:
    input:
        assembly = rules.assembly.output.output_assembly,
        datefile = rules.update_MLST.output.datefile
    output:
        # "%s/{sample}/MLST/{sample}.tsv" %output_folder
        mlst_file = "%s/{sample}/mlst/{assembler}_mlst.tsv" %output_folder,
        done = temp("%s/{sample}/mlst/{assembler}.done" %output_folder)
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

        echo "MLST completed succesfully for {wildcards.sample} with {wildcards.assembler} assembly" > {output.done}
    	"""

# Rule Kleborate:
# Runs Kleborate characterising virulence and resistance in pathogen assemblies
rule kleborate:
    input:
        assembly = rules.assembly.output.output_assembly
    output:
        kleborate_outdir = directory("%s/{sample}/kleborate/{assembler}" %output_folder),
        done = temp("%s/{sample}/kleborate/{assembler}.done" %output_folder)
    params:
        options = lambda wildcards: species_configs[sample_to_organism[wildcards.sample]]["analyses_to_run"]["kleborate"]["options"]
    conda:
        "../envs/kleborate.yaml"
    log:
    	stdout = "Logs/{sample}/Kleborate_{assembler}.log"
    message:
    	"[Kleborate]: Running Kleborate on {wildcards.assembler} assembly from {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.kleborate_outdir})

        cmd="kleborate --assemblies {input.assembly} --outdir {output.kleborate_outdir} {params.options}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "Kleborate completed succesfully for {wildcards.sample} with {wildcards.assembler} assembly" > {output.done}
    	"""

# Rule meningotype
# Runs meningotype on SPAdes assembled contigs to perform serotyping of N. Meningmeningitidis contigs
rule meningotype:
    input:
        assembly = rules.assembly.output.output_assembly
    output:
        meningotype = "%s/{sample}/meningotype/{assembler}_meningotype.tsv" %output_folder,
        done = temp("%s/{sample}/meningotype/{assembler}.done" %output_folder)
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
    
        echo "Meningotype completed succesfully on {wildcards.sample} with {wildcards.assembler} assembly" > {output.done}
    	"""

rule seqsero2:
    input:
        R1 = lambda wildcards: sample_to_illumina[wildcards.sample][0],
        R2 = lambda wildcards: sample_to_illumina[wildcards.sample][1],
    output:
        out_dir = directory("%s/{sample}/seqsero2" %output_folder),
        done = temp("%s/{sample}/seqsero2/seqsero2.done" %output_folder)
    threads:
        max(1, workflow.cores * 0.3333333)
    conda:
        "../envs/seqsero2.yaml"
    log:
        stdout = 'Logs/{sample}/seqsero2.log'
    message:
        "[seqsero2]: Running seqsero2 on {wildcards.sample}"
    shell:
        """
        mkdir -p {output.out_dir}

        cmd="SeqSero2_package.py -m k -t 2 -b mem -i {input.R1} {input.R2} -d {output.out_dir} -n {wildcards.sample} -p {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        echo "seqsero2 completed succesfully for {wildcards.sample}" > {output.done}
        """

rule sistr:
    input:
        assembly = rules.assembly.output.output_assembly,
        serovarlist = rules.fetch_Senterica_Serovar.output.source
    output:
        sistr_tab = "%s/{sample}/sistr/{assembler}_sistr.tab" %output_folder,
        gmlst_profile = "%s/{sample}/sistr/{assembler}_cgmlst_profiles.csv" %output_folder,
        allele_results = "%s/{sample}/sistr/{assembler}_allele-results.json" %output_folder,
        done = temp("%s/{sample}/sistr/{assembler}.done" %output_folder)
    threads:
        max(1, workflow.cores * 0.3333333)
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

        echo "sistr completed succesfully on {wildcards.sample} with {wildcards.assembler} assembly" > {output.done}
        """
