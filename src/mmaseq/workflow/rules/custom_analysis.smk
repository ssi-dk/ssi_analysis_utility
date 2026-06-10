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
        identified_variants = f"{outdir}/{{sample}}/raw/deletion_identifier/deletion_identifier_{{assembler}}_{{database}}.tsv"
    conda:
        ENVS_DIR / "py_utls.yaml"
    log:
        stdout = f"{logdir}/deletion_identifier_{{assembler}}_{{database}}_{{sample}}.log"
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
