
rule setup_PlasmidFinder:
    conda:
        config["analysis_settings"]["plasmidfinder"]["yaml"]
    output: 
        database = directory(f'{database_path}/{config["analysis_settings"]["plasmidfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_PlasmidFinder.log'
    message:
        "[setup_PlasmidFinder]: Setting up PlasmidFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[plasmidfinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the plasmidfinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        date -I > {output.database}/creation.date
        """


rule setup_ResFinder:
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["resfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_ResFinder.log'
    message:
        "[setup_ResFinder]: Setting up ResFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"

            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[resfinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the resfinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        date -I > {output.database}/creation.date
        """


rule setup_PointFinder:
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["pointfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_PointFinder.log'
    message:
        "[setup_PointFinder]: Setting up PointFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1

            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[pointfinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the pointfinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        date -I > {output.database}/creation.date
        """


rule setup_DisinFinder:
    conda:
        config["analysis_settings"]["resfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["disinfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_DisinFinder.log'
    message:
        "[setup_DisinFinder]: Setting up DisinFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/disinfinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[disinfinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the disinfinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        date -I > {output.database}/creation.date
        """


rule setup_VirulenceFinder:
    conda:
        config["analysis_settings"]["virulencefinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["virulencefinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_VirulenceFinder.log'
    message:
        "[setup_VirulenceFinder]: Setting up VirulenceFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[virulencefinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the virulencefinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        
        date -I > {output.database}/creation.date
        """


rule setup_SerotypeFinder:
    conda:
        config["analysis_settings"]["serotypefinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["serotypefinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_SerotypeFinder.log'
    message:
        "[setup_SerotypeFinder]: Setting up SerotypeFinder database"
    shell:
        """
        git clone https://bitbucket.org/genomicepidemiology/serotypefinder_db.git {output.database} > {log.stdout} 2>&1

        for fasta in $(find {output.database} -iname '*.fsa'); do
            idx_prefix={output.database}/$(basename $fasta .fsa)
            cmd="kma index -i $fasta -o $idx_prefix"
            
            echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
            eval $cmd >> {log.stdout} 2>&1
            
            if [ -z $idx_prefix.comb.b ]; then
                echo '[serotypefinder_db]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the serotypefinder_db has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
            fi
        done
        """


rule setup_AMRFinder:
    conda:
        config["analysis_settings"]["amrfinder"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["amrfinder"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_AMRFinder.log'
    message:
        "[setup_AMRFinder]: Setting up AMRFinderPlus database"
    shell:
        """
        cmd="amrfinder_update --database $(dirname {output.database})"
            
        echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        date -I > {output.database}/creation.date
        """

rule setup_EcoliKmerAligner:
    conda:
        config["analysis_settings"]["kmeraligner"]["yaml"]
    params:
        db_prefix = 'ecoligenes'
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["kmeraligner"]["database"]}')
    log:
        stdout = f'Logs/Databases/setup_EcoliKmerAligner.log'
    message:
        "[setup_EcoliKmerAligner]: Setting up EcoliKmerAligner"
    shell:
        """
        mkdir -p {output.database}
        cmd="curl https://raw.githubusercontent.com/ssi-dk/ecoli_fbi/refs/heads/main/db/ecoligenes/ecoligenes.fsa -o {output.database}/{params.db_prefix}.fsa"

        echo "Executing command:\n$cmd\n" > {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        fasta={output.database}/ecoligenes.fsa

        idx_prefix={output.database}/{params.db_prefix}
        cmd="kma index -i $fasta -o $idx_prefix"

        echo "Executing command:\n$cmd\n" >> {log.stdout} 2>&1
        eval $cmd >> {log.stdout} 2>&1

        if [ -z $idx_prefix.comb.b ]; then
            echo '[setup_EcoliKmerAligner]: ERROR - $idx_prefix.comb.b was not created during KMA indexing. This likely means that the EcoliKmerAligner database has changed. Post this message on our Github repository!' 2>&1 >> {log.stdout}
        fi

        date -I > {output.database}/creation.date
        """

rule update_MLST:
    conda:
        config["analysis_settings"]["mlst"]["yaml"]
    output:
        datefile = f'{database_path}/mlst/creation.date'
    log:
        stdout = f'Logs/Databases/update_MLST.log'
    message:
        "[update_MLST]: Updating MLST databases."
    shell:
        """
        DIR=$(which mlst)
        MLSTDIR="$DIR/../../db/pubmlst"


        mlst-download_pub_mlst -d $MLSTDIR  > {log.stdout} 2>&1
        mlst-make_blast_db >> {log.stdout} 2>&1 && date -I > {output.datefile}
        """

rule setup_CdiffToxin:
    conda:
        config["analysis_settings"]["Clostridioides_difficile_db"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["Clostridioides_difficile_db"]["database"]}/Toxin')
    params:
        accession_loci = "AM180355.1:tcdA,tcdB,tcdC AF271719.1:cdtA,cdtB",
        db_toxin = "Cdiff_Toxin"
    log:
        stdout = f'Logs/Databases/setup_CdiffToxin.log'
    message:
        "[setup_CdiffToxin]: Setting up C. difficile toxin database"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.database}
        > {output.database}/{params.db_toxin}.txt
        > {output.database}/{params.db_toxin}.bed6
        > {output.database}/{params.db_toxin}.fasta

        for item in {params.accession_loci}; do
            acc=$(echo $item | cut -d':' -f1)
            loci=$(echo $item | cut -d':' -f2 | tr ',' ' ')

            cmd="python workflow/scripts/genbank_fetcher.py \
                -a $acc --locus $loci \
                -o {output.database}/{params.db_toxin}.txt \
                --bed {output.database}/{params.db_toxin}.bed6 \
                --fasta {output.database}/{params.db_toxin}.fasta \
                --merge 500 --append"

            echo -e \"\nExecuting command:\n$cmd\n\" >> {log.stdout}
            eval $cmd >> {log.stdout} 2>&1
        done

        fasta={output.database}/{params.db_toxin}.fasta

        cmd="samtools faidx $fasta"
        echo "Executing command:\n$cmd\n" >> {log.stdout}
        eval $cmd >> {log.stdout} 2>&1

        idx_prefix={output.database}/{params.db_toxin}
        cmd="kma index -i $fasta -o $idx_prefix"
        echo "Executing command:\n$cmd\n" >> {log.stdout}
        eval $cmd >> {log.stdout} 2>&1

        if [ ! -f "$idx_prefix.comb.b" ]; then
            echo '[CdiffToxin]: ERROR - $idx_prefix.comb.b not created' >> {log.stdout}
        fi
        """

rule setup_CdiffTRST:
    conda:
        config["analysis_settings"]["Clostridioides_difficile_db"]["yaml"]
    output:
        database = directory(f'{database_path}/{config["analysis_settings"]["Clostridioides_difficile_db"]["database"]}/TRST')
    params:
        TR_repeat_sequences = "TR6 TR10",
        TR_repeat_types = "TR6 TR10 TRST"
    log:
        stdout = f'Logs/Databases/setup_CdiffTRST.log'
    message:
        "[setup_CdiffTRST]: Downloading TRST repeat sequences and types"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.database}

        # TR-type repeat sequence downloads
        for TR in {params.TR_repeat_sequences}; do
            cmd="curl -fSL https://raw.githubusercontent.com/ssi-dk/cdiff_fbi/refs/heads/raah_dev/db/TRST/${{TR}}_repeat_sequences.fa \
                -o {output.database}/${{TR}}_repeat_sequences.fa"
            echo "Executing command:\n$cmd\n" >> {log.stdout}
            eval $cmd >> {log.stdout} 2>&1
        done

        # TR-type repeat types downloads
        for TR in {params.TR_repeat_types}; do
            cmd="curl -fSL https://raw.githubusercontent.com/ssi-dk/cdiff_fbi/refs/heads/raah_dev/db/TRST/${{TR}}_repeat_types.txt \
                -o {output.database}/${{TR}}_repeat_types.txt"
            echo "Executing command:\n$cmd\n" >> {log.stdout}
            eval $cmd >> {log.stdout} 2>&1
        done
        """

rule setup_all_databases:
    input:
        rules.setup_EcoliKmerAligner.output.database,
        rules.setup_PlasmidFinder.output.database,
        rules.setup_ResFinder.output.database,
        rules.setup_PointFinder.output.database,
        rules.setup_DisinFinder.output.database,
        rules.setup_VirulenceFinder.output.database,
        rules.setup_AMRFinder.output.database,
        rules.setup_CdiffToxin.output.database,
        rules.setup_CdiffTRST.output.database
