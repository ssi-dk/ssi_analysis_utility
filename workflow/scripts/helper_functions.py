import pandas as pd
import os
import yaml
from typing import Dict, Tuple, Set, List
import warnings
import sys

def determine_sample_configs(samplesheet, config_dir, enable_defaults):
    print("Determining sample configurations")
    # Create a dict for sample names and dict files
    sample_configs = {}

    # Iterate samplesheet and pair samples with configurations
    for sample, cfg in zip(samplesheet["sample_name"], samplesheet["config"]):

        # Deduce configuration file from samplesheet
        cfg_path = f"{config_dir}/{cfg}"

        # Handle missing configuration file
        if not os.path.exists(cfg_path):
            print(f"Warning: Config file specified in samplesheet {cfg} does not exist in {config_dir}!")
            cfg_path = None

            # Check whether to deploy default configurations
            if enable_defaults:
                default_path = f"{config_dir}/default.yaml"

                # Ensure that default file exists and use it
                if os.path.exists(default_path):
                        print("Using default.yaml instead")
                        cfg_path = default_path
                else:
                    print(f"Warning: Default configuration file is missing, please recreate it to enable default analysis: {default_path}")


        # Read sample configrations
        if cfg_path is not None:
            #print(f"Configuration file {cfg_path} found for {sample}") # As log_debug
            with open(cfg_path, "r") as config_file:
                sample_configs[sample] = yaml.safe_load(config_file)

        # Warn user of no configuration is included
        else:
            print(f"No configuration file was specified for sample {sample}.\nSkipping!")
    
    # Ensure that there are indeed sample configurations
    if len(sample_configs) == 0:
        sys.exit("No sample configuration files found. Ensure that the `config` column of the samplesheet is correctly filled.")

    #print(f"Returning sample_configs with length {len(sample_configs)} as {type(sample_configs)}") # As log_debug
    return(sample_configs)


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



def list_results(sample_configs, results_catalogue, output_folder):
    print("Determining results files to create")
    results_files = []

    # Matching each sample-module configurations to the results catalogue
    for sample, modules in sample_configs.items():

        #print(f"Inspecting modules for {sample}") # AS Log_debug
        
        # Extracting module specific configs
        for module in modules:
            configs = modules.get(module)

            # Ensure that there are results file for module
            if module in results_catalogue.keys():
                results_dir = f"{output_folder}/{sample}/{module}"
                results_per_module = results_catalogue.get(module)

                # Ensuring that module results are handles as same object
                if type(results_per_module) is str:
                    results_per_module = [results_per_module]

                # Iterating over individual results files
                for module_results in results_per_module:
                    if type(module_results) is list:
                        if len(module_results) > 1:
                            sys.exit(f"Module results out of bounds: {module_results}")
                        module_results = module_results[0]

                    results_file = f"{results_dir}/{module_results}"

                    # Update results file with configuration
                    if type(configs) is dict:
                        configs["sample"] = sample

                        # Map ALL configs to current result file regardless of being mapable
                        try:
                            # Sloppily map configs to results file    
                            results_file = f"{results_dir}/{module_results.format(**configs)}"
                        except KeyError as e:
                            # Ignore when sloppy mapping fails
                            pass

                    # Collect results files
                    results_files.append(results_file)


    # Ensure that we indeed have results files
    if len(results_files) < 1:
        sys.exit("Warning: No result files detected.")

    return(results_files)
