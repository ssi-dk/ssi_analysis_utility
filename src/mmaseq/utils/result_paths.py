#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from collections import defaultdict

# Import from utils
from .helpers import deconvolute_path, read_results_catalogue


def define_module_results_file(outdir, sample, module, results_catalogue, sample_configs, raw):
    """
    Defines the expected result file paths for a specific module and sample.

    Args:
        outdir (str): Output directory path.
        sample (str): Sample name.
        module (str): Module name.
        results_catalogue (dict): Dictionary containing result file templates for each module.
        sample_configs (dict): Dictionary containing configurations for each sample.

    Returns:
        list: List of Path objects representing the expected result file paths.
    """
    # Define prefix for result file paths
    prefix = Path(f"{outdir}/{sample}/{module}").expanduser()
    if raw:
        prefix = Path(f"{outdir}/{sample}/raw/{module}").expanduser()
    

    # Define and normalise expected result file names
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


def define_all_result_files(outdir, sample_configs, results_catalogue, raw):
    """
    Defines all expected result file paths for all samples and modules.

    Args:
        outdir (str): Output directory path.
        sample_configs (dict): Dictionary containing configurations for each sample.
        results_catalogue (dict): Dictionary containing result file templates for each module.

    Returns:
        defaultdict: Nested dictionary with sample -> module -> list of Path objects.
    """
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
            if not isinstance(result_strings, list):
                result_strings = [result_strings]

            # Extract the configurations as keywords
            configs = modules.get(mod)

            # Streamline object as dict
            if not isinstance(configs, dict):
                configs = dict()  # Not a dict means no keywords

            result_files = define_module_results_file(outdir, sample, mod, results_catalogue, sample_configs, raw)

            all_result_files[sample].update({mod: result_files})

    return all_result_files


def determine_rule_output(outdir, sample, module, results_catalogue, sample_configs):
    """
    Determines the output paths for a rule based on the results catalogue and configurations.

    Args:
        outdir (str): Output directory path.
        sample (str): Sample name.
        module (str): Module name.
        results_catalogue (dict): Dictionary containing result file templates for each module.
        sample_configs (dict): Dictionary containing configurations for each sample.

    Returns:
        list: List of output path strings for the rule.
    """
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

    return [f"{outdir}/{{sample}}/{{module}}/{path}" for path in module_result_path]


