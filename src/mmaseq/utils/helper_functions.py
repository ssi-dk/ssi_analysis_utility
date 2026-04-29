import pandas as pd
import os
import yaml
from typing import Dict, Tuple, Set, List
import warnings
import sys
from pathlib import Path
from collections import defaultdict
import itertools

def read_results_catalogue(results_catalogue_path):
    with open(results_catalogue_path, "r") as catalogue_file:
        try:
            results_catalogue = yaml.safe_load(catalogue_file)
        except:
            print("read_results_catalogue(results_catalogue_path): LAZY DEVELOPPERS... Fill out except!!!")
            results_catalogue = yaml.safe_load(results_catalogue_path)
    
    return(results_catalogue)


def determine_sample_configs(samplesheet, config_dir):
    # Create a dict for sample names and dict files
    sample_configs = {}

    # Iterate samplesheet and pair samples with configurations
    for sample, cfg in zip(samplesheet.index, samplesheet["config"]):

        # Deduce configuration file from samplesheet
        cfg_path = f"{config_dir}/{cfg}"

        # Handle missing configuration file
        if not os.path.isfile(cfg_path):
            print(f"Warning: Config file specified in samplesheet {cfg} does not exist in {config_dir}!")
            cfg_path = None

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


def deconvolute_path(template, configs):
    if not isinstance(configs, dict):
        return [template]

    # Convert multiple arguments into exhaustive parallel lists
    options = list(configs.keys())
    arguments = [arg if isinstance(arg, (list, tuple)) else [arg] for arg in (configs[k] for k in options)]

    exhaustive_templates = list()
    # Iteratively generate exhaustive dicts for handling multiple option/arugment relationships
    for combo in itertools.product(*arguments):
        paired = dict(zip(options, combo))
        # Laizily map option/arguments where applicatble into expected result file names
        if len(configs) > 0:
            fname = template.format(**paired)
            exhaustive_templates.append(fname)
                
        # Handle missing options
        else:
            exhaustive_templates.append(template)

    return exhaustive_templates


def define_module_results_file(outdir, sample, module, results_catalogue, sample_configs):

    prefix = Path(f"{outdir}/{sample}/{module}").expanduser()

    module_result_files = list()

    # Define and normalise expected reult file names
    result_strings = results_catalogue.get(module)
    if not isinstance(result_strings, (list, tuple)):
        result_strings = [result_strings]

    # Define module configurations
    sample_cfg = sample_configs.get(sample)
    configs = sample_cfg.get(module)

    # Define container for results
    module_result_files = []

    # Iterate over multiple expected output results files
    for template in result_strings:
        deconvoluted_paths = deconvolute_path(template, configs)

        for path in deconvoluted_paths:
            module_result_files.append(prefix / path)

    return module_result_files



def define_all_result_files(outdir, sample_configs, results_catalogue):

    # Define carrier object
    all_result_files = defaultdict(dict)

    # Iterate over individual sample configurations
    for sample, modules in sample_configs.items():
        
        # Iterate over individual modules
        for mod in modules.keys():

            # Ensure that module exists in results catalogue
            if mod not in results_catalogue.keys():
                continue
            
            # Define Result file strings
            result_strings = results_catalogue.get(mod)

            # Streamline object as list
            if type(result_strings) is not list:
                result_strings = [result_strings]

            # Extract the configurations as keywords
            configs = modules.get(mod)

            # Streamline object as dict
            if type(configs) is not dict:
                configs = dict() # Not a dict means no keywords

            result_files = define_module_results_file(outdir, sample, mod, results_catalogue, sample_configs)

            all_result_files[sample].update({mod: result_files})

    return all_result_files


def unpivot_results(sample, module, file, results):

    # Convert index to row number column starting from row 1
    results.index += 1
    results = results.reset_index(names = "Row")

    # Generate long list format of results file
    results_long = results.melt(
        id_vars = "Row",
        var_name = "Column",
        value_name = "Value"
        )

    # add columns in a single assignment (faster than multiple insert calls)
    results_long[["Sample", "Module", "File"]] = [sample, module, file.name]
    
    # if you want a specific column order:
    cols = ["Sample", "Module", "File", "Row", "Column", "Value"]
    results_long = results_long[cols]

    return results_long


def generate_long_results(all_result_files):
    all_sample_results = list()

    for sample, modules in all_result_files.items():

        for mod, files in modules.items():

            long_mod_results = list()

            for file in files:

                try:
                    sample_results = pd.read_csv(file, sep = "\t", index_col = False)
                except pd.errors.EmptyDataError as e:
                    print(f"Results file {file.name} for {sample} is empty. Skipping!")
                    continue

                # Determine whether the long table format is allready observed
                sample_long = sample_results
                if not {"Sample", "Module", "File", "Row", "Column", "Value"}.issubset(sample_results.columns):
                    sample_long = unpivot_results(sample, mod, file, sample_results)

                all_sample_results.append(sample_long)

    return pd.concat(all_sample_results, ignore_index = True)


def determine_rule_output(outdir, sample, module, results_catalogue, sample_configs):

    result_strings = results_catalogue.get(module)
    if not isinstance(result_strings, (list, tuple)):
        result_strings = [result_strings]

    # Define module configurations
    sample_cfg = sample_configs.get(sample)
    configs = sample_cfg.get(module)

    # Define container for results
    module_result_path = []

    # Iterate over multiple expected output results files
    for template in result_strings:
        deconvoluted_paths = deconvolute_path(template, configs)

        for path in deconvoluted_paths:
            module_result_path.append(path)

    return ["%s/{sample}/{module}/" %outdir + path for path in module_result_path]