rule spades:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = "%s/{sample}/spades/{sample}.fasta" %output_folder,
        tool_version = "%s/{sample}/spades/spades_version.txt" %output_folder,
    conda:
        "../envs/shovill.yaml"
    log:
        stdout = "Logs/Assemblies/{sample}_spades.log"
    threads:
        max(1, workflow.cores * 0.66667)
    priority: 2
    message:
        "[SPAdes]: Assemblying {wildcards.sample} using SPAdes with {threads} thread(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        outdir=$(dirname {output.assembly})
        cmd="spades.py -1 {input.R1} -2 {input.R2} --threads {threads} --isolate -o $outdir"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="mv $outdir/contigs.fasta {output.assembly}"
        echo "### SPAdes Done!###\nExecuting command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # 2) extract tool version
        version_cmd="spades.py -v"
        date_cmd="date -I"
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}
        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        # Write "<accessions>\t<date>" to the version file
        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
        """


rule skesa:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = "%s/{sample}/skesa/{sample}.fasta" %output_folder,
        tool_version = "%s/{sample}/skesa/skesa_version.txt" %output_folder,
    conda:
        "../envs/shovill.yaml"
    log:
        stdout = "Logs/Assemblies/{sample}_Skesa.log"
    threads:
        max(1, workflow.cores * 0.66667)
    priority: 2
    message:
        "[Skesa]: Assemblying {wildcards.sample} using Skesa with {threads} core(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        cmd="skesa --reads {input.R1},{input.R2} --contigs_out {output.assembly} --cores {threads}"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # 2) extract tool version
        version_cmd="skesa -v"
        date_cmd="date -I"
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}
        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        # Write "<accessions>\t<date>" to the version file
        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
        """


rule shovill:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"]
    output:
        assembly = "%s/{sample}/shovill/{sample}.fasta" %output_folder,
        tool_version = "%s/{sample}/shovill/shovill_version.txt" %output_folder,
    conda:
        "../envs/shovill.yaml"
    log:
        stdout = "Logs/Assemblies/{sample}_Shovill.log"
    threads:
        max(1, workflow.cores * 0.66667)
    priority: 2
    message:
        "[Shovill]: Assemblying {wildcards.sample} using Shovill with {threads} CPU(s). This may take some time!\nInspect {log.stdout} for more details!"
    shell:
        """
        mkdir -p $(dirname {output.assembly})
        outdir=$(dirname {output.assembly})

        cmd="shovill --R1 {input.R1} --R2 {input.R2} --outdir $outdir/ --force --cpus {threads}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        cmd="mv $outdir/contigs.fa {output.assembly}"
        echo "### Shovil Done!###\nExecuting command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # 2) extract tool version
        version_cmd="shovill --version"
        date_cmd="date -I"
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}
        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        # Write "<accessions>\t<date>" to the version file
        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
        """


rule assembly:
    input:
        input_assembly = "%s/{sample}/{assembler}/{sample}.fasta" %output_folder,
        input_assembly_version = "%s/{sample}/{assembler}/{assembler}_version.txt" %output_folder,
    output:
        output_assembly = "%s/{sample}/Assemblies/{sample}_{assembler}.fasta" %output_folder,
        output_assembly_version = "%s/{sample}/Assemblies/{sample}_{assembler}_version.txt" %output_folder,
    log:
        stdout = "Logs/Assemblies/{sample}_{assembler}_assembly.log"
    shell:
        """
        mkdir -p $(dirname {output.output_assembly})

        cmd="cp {input.input_assembly} {output.output_assembly}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1 
        
        cmd_version="cp {input.input_assembly_version} {output.output_assembly_version}"

        echo "Executing command:\n$cmd_version\n" >> {log.stdout} 2>&1
        eval $cmd_version >> {log.stdout} 2>&1
        """
