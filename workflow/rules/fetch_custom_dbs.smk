rule fetch_genbank:
  output:
    fasta = temp("%s/GenBank/{accession}-{loci}.fasta" %database_path),
    bed = temp("%s/GenBank/{accession}-{loci}.bed6" %database_path),
    meta = temp("%s/GenBank/{accession}-{loci}.txt" %database_path)
  conda:
    config["analysis_settings"]["Clostridioides_difficile_db"]["yaml"]
  log:
    stdout = 'Logs/Databases/fetch_genbank_{accession}-{loci}.log'
  message:
    "[fetch_genbank]: Fetching {wildcards.accession} from Genbank"
  shell:
    """
    cmd="python workflow/scripts/genbank_fetcher.py -a {wildcards.accession} --locus $(echo {wildcards.loci} | tr _ ' ') -o {output.meta} --bed {output.bed} --fasta {output.fasta} --merge 500 --append"

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


rule fetch_cdifftoxins:
  input:
    fastas = ["%s/GenBank/AM180355.1-tcdA_tcdB_tcdC.fasta" %database_path, "%s/GenBank/AF271719.1-cdtA_cdtB.fasta" %database_path],
    beds = ["%s/GenBank/AM180355.1-tcdA_tcdB_tcdC.bed6" %database_path, "%s/GenBank/AF271719.1-cdtA_cdtB.bed6" %database_path]
  output:
    fasta = "%s/custom/Cdiff_Toxins.fasta" %database_path,
    bed = "%s/custom/Cdiff_Toxins.bed6" %database_path,
  shell:
    """
    cat {input.fastas} > {output.fasta}

    header=$(grep "#" -h {input.beds} | uniq)
    genes=$(grep "#" -vh {input.beds} | uniq)
    echo "$header\n$genes" > {output.bed}
    """


rule fetch_chtyper_db:
  output:
    source = "%s/custom/fumCH_db.fasta" %database_path
  conda:
    "../envs/fetch.yaml"
  log:
    stdout = "Logs/Databases/setup_chtyper_database.log"
  message:
    "[fetch_chtyper_db]: Downloading custom database for CHtyper"
  shell:
    """
    outdir=$(dirname {output.source})
    mkdir -p $outdir
    
    cmd="curl https://bitbucket.org/genomicepidemiology/chtyper_db/raw/654ca48d250e0a69c6c06b4be5a96d807b23f806/fimH.fsa -o $outdir/fimH.fsa"

    echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
    eval $cmd >> {log.stdout} 2>&1
    
    cmd="curl https://bitbucket.org/genomicepidemiology/chtyper_db/raw/654ca48d250e0a69c6c06b4be5a96d807b23f806/fumC.fsa -o $outdir/fumC.fsa"
    eval $cmd >> {log.stdout} 2>&1 

    cmd="cat $outdir/fimH.fsa $outdir/fumC.fsa > {output.source}"
    eval $cmd >> {log.stdout} 2>&1
    """
