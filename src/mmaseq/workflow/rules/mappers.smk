rule custom_kmeralignment:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_custom_kmeraligner_index.output.names
    params:
        prefix_out = "%s/{sample}/kmeraligner/{database}" %outdir,
        prefix_db = rules.setup_custom_kmeraligner_index.params.prefix    
    output:
        results = "%s/{sample}/kmeraligner/{database}.res" %outdir,
        sam = temp("%s/{sample}/samtools/{database}.sam" %outdir),
        matrix = temp("%s/{sample}/kmeraligner/{database}.mat.gz" %outdir),
        tool_version = "%s/{sample}/kmeraligner/{database}_kmaalign_version.txt" %outdir,
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = "%s/{sample}/custom_kmeralignment_{database}.log" %logdir
    message:
        "[kmeraligner]: Running KMA for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.results})
        mkdir -p $(dirname {output.sam})

        cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix_out} -t_db {params.prefix_db} -na -nc -nf -sam 4 -matrix > {output.sam}"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # 2) create version file with date
        version_cmd="kma -v"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
        """

rule custom_kmerconsensus:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_custom_kmeraligner_index.output.names
    params:
        prefix_out = "%s/{sample}/kmerconsensus/{database}" %outdir,
        prefix_db = rules.setup_custom_kmeraligner_index.params.prefix,
    output:
        results = temp("%s/{sample}/kmerconsensus/{database}.res" %outdir),
        seq = "%s/{sample}/kmerconsensus/{database}.fsa" %outdir,
        aln = temp("%s/{sample}/kmerconsensus/{database}.aln" %outdir),
        tool_version = "%s/{sample}/kmerconsensus/{database}_kmaconsensus_version.txt" %outdir,
    conda:
        ENVS_DIR / "kmeraligner.yaml"
    log:
        stdout = "%s/{sample}/custom_kmerconsensus_{database}.log" %logdir
    message:
        "[kmerconsensus]: Running KMA for {wildcards.database} on {wildcards.sample}"
    shell:
        """
        mkdir -p $(dirname {output.seq})

        cmd="kma -ipe {input.R1} {input.R2} -o {params.prefix_out} -t_db {params.prefix_db} -nf -ref_fsa"
        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        # 2) create version file with date
        version_cmd="kma -v"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
        """

rule custom_bowtie2alignment:
    input:
        R1 = lambda wc: samplesheet.loc[wc.sample, "read1"],
        R2 = lambda wc: samplesheet.loc[wc.sample, "read2"],
        database = rules.setup_custom_bowtie2_index.output.bt2_1 # just locate one of the bt2 files to activate the db_setup
    params:
       options = lambda wc: sample_configs[wc.sample]["custom_bowtie2alignment"]["options"]
    output:
        sam = temp("%s/{sample}/bowtie2/{database}.sam" %outdir)
    threads: workflow.cores - 1 - (workflow.cores - 1) % 2
    priority: 1
    conda:
        ENVS_DIR / "bowtie2.yaml"
    log:
        stdout = "%s/{sample}/custom_bowtie2_{database}.log" %logdir
    message:
        "[bowtie2aligner]: Running Bowtie2 for {wildcards.database} on {wildcards.sample} using {threads} thread(s)"
    shell:
        """
        mkdir -p $(dirname {output.sam})

        db_prefix="{input.database}"
        db_prefix="${{db_prefix%.1.bt2}}"
        
        cmd="bowtie2 -1 {input.R1} -2 {input.R2} -q -S {output.sam} {params.options} -x $db_prefix --threads {threads}"
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
        results = "%s/{sample}/custom_blaster/blast_{assembler}_{database}.tsv" %outdir,
        tool_version = "%s/{sample}/custom_blaster/blast_{assembler}_{database}_version.txt" %outdir,
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
    
        # 2) create version file with date
        version_cmd=" blastn -version|head -1"
        date_cmd="date -I"
            
        echo -e "Executing command:\n$version_cmd\n$date_cmd\n" >> {log.stdout}

        version_str="$(eval "$version_cmd" 2>> {log.stdout})"
        date_str="$(eval "$date_cmd" 2>> {log.stdout})"

        printf '%s\t%s\n' "$version_str" "$date_str" > {output.tool_version}
        """
        
rule assembly_minimap2:
    input:
        assembly = rules.assembly.output.output_assembly,      # {sample}, {assembler}
        database = rules.fetch_genbank.output.fasta            # {sample}, {database}
    params:
        options = lambda wc: sample_configs[wc.sample]["assembly_minimap2"]["options"]
    output:
        results = temp(f"{outdir}/{{sample}}/minimap2/{{assembler}}_{{database}}.sam")
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