import pandas as pd
import os
import yaml
from typing import Dict, Tuple, Set, List
import warnings
import sys
import itertools
from copy import deepcopy


def read_results_catalogue(results_catalogue_path):
    with open(results_catalogue_path, "r") as catalgoue_file:
        try:
            results_catalogue = yaml.safe_load(catalgoue_file)
        except:
            print("read_results_catalogue(results_catalogue_path): LAZY DEVELOPPERS... Fill out except!!!")
            yaml.safe_load(results_catalogue_path)

    return(results_catalogue)


def determine_sample_configs(samplesheet, config_dir, enable_defaults):
    #print("Determining sample configurations") #DEBUG msg
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
) -> Tuple[Dict[str, str], Dict[str, List[str]], Dict[str, str], Dict[str, str]]:
    """
    read_base and assembly_base are the paths defined in the config files for the input data

    read_path: examples/Dataset/reads
    assembly_path: examples/Dataset/assembly_path

    with the function returning 
    #{'ERR3528110': ['examples/Dataset/reads/ERR3528110_1.fastq.gz', 'examples/Dataset/reads/ERR3528110_2.fastq.gz'],
    """

    samples = {}
    sample_to_illumina = {} 
    sample_to_nanopore = {}
    sample_to_assembly = {}

    for idx, row in samplesheet.iterrows():
        sample = row[sample_col] #the key -> 'ERR3528110'

        samples[sample] = sample 
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

    return samples, sample_to_illumina, sample_to_nanopore, sample_to_assembly


def map_configs_to_results(configs, results_dir, results_file):
    """
    Expand list-valued configuration fields into all possible combinations
    and format the result file pattern for each concrete configuration set.

    Example:
        configs = { "assemblers": ["spades","skesa"], "sample": "ERR142064" }
        results_file = "{assemblers}_repeat_types.tsv"

    Produces:
        [
            ".../spades_repeat_types.tsv",
            ".../skesa_repeat_types.tsv"
        ]
    """
    # Work on a copy to avoid mutating the original dictionary
    cfg = deepcopy(configs)

    # Prepare to expand every config key; scalar values become 1-item lists
    keys = list(cfg.keys())
    value_lists = []
    for k in keys:
        v = cfg[k]
        if isinstance(v, list):
            # Already a list → keep as is
            value_lists.append(v)
        else:
            # Scalar → wrap into list so itertools.product works uniformly
            value_lists.append([v])

    results = []

    # Cartesian product over all config values
    # → yields one fully concrete mapping per file to generate
    for combo in itertools.product(*value_lists):
        mapping = dict(zip(keys, combo))  # map keys to this particular combination

        try:
            # Attempt to apply the mapping to the filename format string
            path = f"{results_dir}/{results_file.format(**mapping)}"
        except KeyError:
            # If the result pattern references a config key that isn't present, silently skip.
            continue

        results.append(path)

    return results


def list_results(sample_configs, results_catalogue, output_folder):
    """
    Construct a full list of output file paths for each sample by:
      - iterating over all sample modules
      - matching module names to the results catalogue
      - expanding list-valued configuration fields (assemblers, repeats, etc.)
      - formatting each results entry using module configs
    """

    results_files = []
    seen = set()  # Used to deduplicate while preserving ordering

    # Iterate over all modules for individual samples
    for sample, modules in sample_configs.items():

        # Iterate over all configurations for any sample specific module
        for sample_module, configs in modules.items():

            # Skip modules that are NOT in the results catalogue
            if sample_module not in results_catalogue:
                continue

            # Construct the directory where results for this module live
            results_dir = f"{output_folder}/{sample}/{sample_module}"

            # Normalise results definitions into a list for uniform iteration
            module_results = results_catalogue.get(sample_module)
            if isinstance(module_results, str):
                module_results = [module_results]

            # Case 1: Current module has configurations that must be mapped to the results file
            if isinstance(configs, dict):

                # Copy configs and add the sample name for sloppy mapping
                configs_copy = deepcopy(configs)
                configs_copy["sample"] = sample

                # Handle multiple module result files
                for results_file in module_results:

                    # Map configuration parameters to results files
                    concrete_paths = map_configs_to_results(
                        configs_copy,
                        results_dir,
                        results_file
                    )

                    # Prevent potential duplications, then record final paths
                    for path in concrete_paths:
                        if path not in seen:
                            seen.add(path)
                            results_files.append(path)

            # Case 2: Current module has *no* configs at all (rare but allowed)
            else:
                # Use raw file patterns with no mapping
                for results_file in module_results:
                    path = f"{results_dir}/{results_file}"
                    if path not in seen:
                        seen.add(path)
                        results_files.append(path)

    # Fail fast if no files were generated at all
    if not results_files:
        sys.exit("Warning: No result files detected.")

    return results_files


def read_tsv(path):
    raw = pd.read_csv(path, sep = "\t")
    file = os.path.basename(path)
    raw["file"] = file
