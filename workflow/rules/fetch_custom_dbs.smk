rule fetch_genbank:
  input:
    meta_file = "%s/../examples/Metadata/{database}_genbank_metafile.tsv" %database_path
  output:
    fasta = "%s/custom/{database}.fasta" % database_path,
    bed = "%s/custom/{database}.bed6" % database_path,
    records = "%s/custom/{database}.txt" % database_path
  params:
    merge = 500
  conda:
    "../envs/DatabaseFetch.yaml"
  log:
    stdout = 'Logs/Databases/fetch_genbank_{database}.log'
  message:
    "[fetch_genbank]: Fetching {wildcards.database} from Genbank"
  shell:
    """
    mkdir -p $(dirname {output.fasta})

    cmd="python workflow/scripts/genbank_fetcher.py --metafile {input.meta_file} --bed {output.bed} --records {output.records} --fasta {output.fasta} --merge {params.merge} --append"

    echo "Executing command:\n$cmd\n" > {log.stdout}
    eval $cmd >> {log.stdout} 2>&1
    """

rule fetch_type_repeat_sequence:
  output:
    seq = "%s/custom/type_repeats/{TR}.fasta" %database_path
  conda:
    "../envs/fetch.yaml"
  log:
    stdout = "Logs/Databases/fetch_type_repeat_sequences_{TR}.log"
  message:
    "[fetch_type_repeat_sequences]: Downloading Type Repeat Sequence Type sequences"
  shell:
    """
    mkdir -p $(dirname {output.seq})

    cmd="curl -fSL https://raw.githubusercontent.com/ssi-dk/cdiff_fbi/refs/heads/raah_dev/db/TRST/{wildcards.TR}_repeat_sequences.fa -o {output.seq}"

    echo "Executing command:\n$cmd\n" > {log.stdout}
    eval $cmd >> {log.stdout} 2>&1
    """


rule fetch_type_repeat_metadata:
  output:
    meta = "%s/custom/{TR}.txt" %database_path
  conda:
    "../envs/fetch.yaml"
  log:
    stdout = "Logs/Databases/fetch_type_repeat_metadata_{TR}.log"
  message:
    "[fetch_type_repeat_metadata]: Downloading Type Repeat Sequence Type metadata"
  shell:
    """
    mkdir -p $(dirname {output.meta})

    cmd="curl -fSL https://raw.githubusercontent.com/ssi-dk/cdiff_fbi/refs/heads/raah_dev/db/TRST/{wildcards.TR}_repeat_types.txt -o {output.meta}"
    
    echo "Executing command:\n$cmd\n" > {log.stdout}
    eval $cmd >> {log.stdout} 2>&1
    """


rule fetch_ecoligenes:
  output:
    source = "%s/custom/ecoligenes.fasta" %database_path
  conda:
    "../envs/fetch.yaml"
  log:
    stdout = "Logs/Databases/setup_ecoligenes_ecoligenes.log"
  message:
    "[fetch_ecoligenes]: Downloading custom database ecoligenes"
  shell:
    """
    mkdir -p $(dirname {output.source})
    cmd="curl https://raw.githubusercontent.com/ssi-dk/ecoli_fbi/refs/heads/main/db/ecoligenes/ecoligenes.fsa -o {output.source}"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    """