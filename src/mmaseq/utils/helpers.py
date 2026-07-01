#!/usr/bin/env python3

import yaml
import itertools


def read_results_catalogue(results_catalogue_path):
    """
    Reads the results catalogue from a YAML file.

    Args:
        results_catalogue_path (str): Path to the YAML file containing the results catalogue.

    Returns:
        dict: The loaded results catalogue dictionary.
    """
    # Read results catalogue from yaml file
    with open(results_catalogue_path, "r") as catalogue_file:
        try:
            results_catalogue = yaml.safe_load(catalogue_file)
        except:
            print("read_results_catalogue(results_catalogue_path): LAZY DEVELOPPERS... Fill out except!!!")
            results_catalogue = yaml.safe_load(results_catalogue_path)
    
    return results_catalogue


def deconvolute_path(template, configs):
    """
    Deconvolutes a path template by substituting configurations into it, handling multiple values.

    Args:
        template (str): Path template string with placeholders.
        configs (dict): Dictionary of configuration options and their values.

    Returns:
        list: List of deconvoluted path strings.
    """
    if not isinstance(configs, dict):
        return [template]

    # Convert multiple arguments into exhaustive parallel lists
    options = list(configs.keys())
    arguments = [arg if isinstance(arg, (list, tuple)) else [arg] for arg in (configs[k] for k in options)]

    exhaustive_templates = list()
    # Iteratively generate exhaustive dicts for handling multiple option/argument relationships
    for combo in itertools.product(*arguments):
        paired = dict(zip(options, combo))
        # Lazily map option/arguments where applicable into expected result file names
        if len(configs) > 0:
            fname = template.format(**paired)
            exhaustive_templates.append(fname)
                
        # Handle missing options
        else:
            exhaustive_templates.append(template)

    return exhaustive_templates