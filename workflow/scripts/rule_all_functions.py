import pandas as pd
import os
import yaml
from typing import Dict, Tuple, Set, List
import warnings

def build_config_lookups(samplesheet, config_dir, sample_col="sample_name", config_col="config"):
    """
    Returns two plain dict lookups:

      options_lookup[sample][analysis] -> options string
      tool_lookup[sample][tool][subtool] -> options string

    """
    options_lookup = {}   # { sample: { analysis: "options" } }
    tool_lookup    = {}   # { sample: { tool: { subtool: "options" } } }

    for _, row in samplesheet.iterrows():
        sample   = str(row[sample_col])
        config_name = str(row[config_col])
        config_path = config_name if os.path.isabs(config_name) else os.path.join(config_dir, config_name)

        # ensure nested dicts exist
        if sample not in options_lookup:
            options_lookup[sample] = {}
        if sample not in tool_lookup:
            tool_lookup[sample] = {}

        config = {}
        if not os.path.exists(config_path):
            warnings.warn(f"Config file not found for sample {sample}: {config_path}. Using empty.")
        else:
            try:
                with open(config_path, "r") as fh:
                    config = yaml.safe_load(fh) or {}
            except Exception as e:
                warnings.warn(f"Could not read YAML for {sample} ({config_path}): {e}. Using empty.")

        # --- analyses_to_run: store only 'options' values ---
        if "analyses_to_run" in config and isinstance(config["analyses_to_run"], dict):
            for analysis_key, settings in config["analyses_to_run"].items():
                if isinstance(settings, dict) and "options" in settings:
                    options_lookup[sample][analysis_key] = settings["options"]

        # --- tool_settings: store only 'options' values ---
        if "tool_settings" in config and isinstance(config["tool_settings"], dict):
            for tool, subtools in config["tool_settings"].items():           # e.g., samtools
                if not isinstance(subtools, dict):
                    continue
                if tool not in tool_lookup[sample]:
                    tool_lookup[sample][tool] = {}
                for subtool, params in subtools.items():                  # e.g., view / sort
                    if isinstance(params, dict) and "options" in params:
                        tool_lookup[sample][tool][subtool] = params["options"]

    return options_lookup, tool_lookup

def sample_read_map(
    samplesheet: pd.DataFrame,
    sample_col: str = "sample_name",
    read_col: str = "Illumina_read_files",
    nanopore_col: str = "Nanopore_read_file",
    assembly_col: str = "assembly_file",
    illumina_read_base: str = "",
    nanopore_read_base: str = "",
    assembly_base: str = "",
) -> Tuple[Dict[str, List[str]], Dict[str, str], Dict[str, str]]:
    """
    read_base and assembly_base are the paths defined in the config files for the input data

    read_path: examples/Dataset/reads
    assembly_path: examples/Dataset/assembly_path

    with the function returning 
    #{'ERR3528110': ['examples/Dataset/reads/ERR3528110_1.fastq.gz', 'examples/Dataset/reads/ERR3528110_2.fastq.gz'],
    """

    sample_to_illumina = {} 
    sample_to_nanopore = {}
    sample_to_assembly = {}

    for idx, row in samplesheet.iterrows():
        sample = row[sample_col] #the key -> 'ERR3528110'

        # --- Illumina ---
        illumina_reads = str(row[read_col]).split(',')
        illumina_full = [os.path.join(illumina_read_base, read) for read in illumina_reads]
        sample_to_illumina[sample] = illumina_full

        # --- Nanopore ---
        Nanopore_reads = str(row[nanopore_col])
        Nanopore_full = os.path.join(nanopore_read_base, Nanopore_reads)
        sample_to_nanopore[sample] = Nanopore_full

        # --- Illumina ---
        Assembly_seq = str(row[assembly_col])
        Assembly_full = os.path.join(assembly_base, Assembly_seq)
        sample_to_assembly[sample] = Assembly_full

    return sample_to_illumina, sample_to_nanopore, sample_to_assembly

def list_results(samplesheet, output_folder, species_config_path):
    results = set()

    for _, row in samplesheet.iterrows():
        sample = str(row["sample_name"])
        config_file = str(row["config"])

        #combine the path from the config.yaml file with information from samplesheet like e.coli.yaml
        species_config_file = os.path.join(species_config_path,config_file)

        if not os.path.exists(species_config_file):
            warnings.warn(f"Config file not found for sample {sample}: {species_config_file}. Skipping.")
            continue

        # Load config to determine analyses to run
        with open(species_config_file, "r") as f:
            try:
                config_data = yaml.safe_load(f)
            except Exception as e:
                warnings.warn(f"Could not read YAML for {sample} ({species_config_file}): {e}")
                config_data = {}
               
        # Require analyses_to_run section
        if "analyses_to_run" not in config_data:
            warnings.warn(f"No 'analyses_to_run' section in {species_config_file}. Skipping {sample}.")
            continue

        analyses = config_data["analyses_to_run"]

        for analysis_name, settings in analyses.items():
            # Case 1: Simple true/false flag - like plasmid finder where no database or other wildcards are currently needed
            if not isinstance(settings, dict):
                if settings:
                    results.add(f"{output_folder}/{sample}/{analysis_name}/{analysis_name}.done")
                continue

            # Case 2: Dict with assemblers and/or databases
            if "assemblers" in settings:
                assemblers = settings["assemblers"]
                # ensure assemblers are a list of strings
                if isinstance(assemblers, str):
                    assemblers = [assemblers]
            else:
                assemblers = None

            if "database" in settings:
                databases = settings["database"]
                # database case (can be single string or list of strings)
                if isinstance(databases, str):
                    databases = [databases]
            else:
                databases = None
            
            if databases is None:
                 # no database case
                if assemblers:
                    for tool in assemblers:
                        results.add(f"{output_folder}/{sample}/{analysis_name}/{tool}.done")
                else:
                    results.add(f"{output_folder}/{sample}/{analysis_name}/{analysis_name}.done")
            else:
                if assemblers:
                    for tool in assemblers:
                        for database in databases:
                            results.add(f"{output_folder}/{sample}/{analysis_name}/{tool}_{database}.done")
                else:
                    for database in databases:
                        results.add(f"{output_folder}/{sample}/{analysis_name}/{database}.done")
    print(sorted(results))

    """
    ['examples/Results/ERR142064/amrfinder/spades.done', 
    'examples/Results/ERR142064/cdiff_repeat_identifier/skesa.done', 
    'examples/Results/ERR142064/cdiff_repeat_identifier/spades.done', 
    .
    .
    .
    'examples/Results/SRR4046826/resfinder/resfinder.done', 
    'examples/Results/SRR4046826/serotypefinder/serotypefinder.done', 
    'examples/Results/SRR4046826/virulencefinder/virulencefinder.done']
    """
    return sorted(results)
