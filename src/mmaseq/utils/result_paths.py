#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from itertools import product
from collections import defaultdict

# Import from utils
from .helpers import deconvolute_path, read_results_catalogue


def define_module_results_file(outdir, sample, module, module_options, raw):
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
    
    # Extract configurations
    read_string = ""
    if module_options.get("reads"):
        read_type = module_options.get("read_type")

        if read_type is not None:
            read_string = f"{read_type}/"
        else:
            print(f"ERROR: Read type not identified for {sample}. Check read paths for your samplesheet!")
    

    databases = [None]
    if "database" in module_options.keys():
        databases = module_options.get("database")
        if not isinstance(databases, list):
            databases = [databases]

    assemblers = [None]
    if "assembler" in module_options.keys():
        assemblers = module_options.get("assembler")

        if not isinstance(assemblers, list):
            assemblers = [assemblers]

    paths = []

    for db, asm in product(databases, assemblers):

        suffix = ""

        if db:
            suffix += f"_{db}"

        if asm:
            suffix += f"_{asm}"

        path = Path(
            f"{outdir}/{sample}/{read_string}{module}/{module}{suffix}.tsv"
        )
        if raw:
            path = f"{outdir}/{sample}/raw/{read_string}{module}/{module}{suffix}.tsv"

        paths.append(path)

    return paths


def define_all_result_files(outdir, sample_configs, raw):
    """
    Defines all expected result file paths for all samples and modules.

    Args:
        outdir (str): Output directory path.
        sample_configs (dict): Dictionary containing configurations for each sample.
        raw (bool): Flag whether to store in raw folder or not
    Returns:
        defaultdict: Nested dictionary with sample -> module -> list of Path objects.
    """
    # Define carrier object
    all_result_files = defaultdict(dict)

    # Iterate over individual sample configurations
    for sample, modules in sample_configs.items():
        
        # Iterate over individual modules
        for mod, opts in modules.items():

            # Ensure that tools are skipped 
            if "ignore" in opts.keys():
                continue

            # Streamline object as dict
            if not isinstance(opts, dict):
                opts = dict()  # Not a dict means no keywords

            result_files = define_module_results_file(outdir, sample, mod, opts, raw)

            all_result_files[sample].update({mod: result_files})

    return all_result_files