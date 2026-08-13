#!/usr/bin/env python3

import pandas as pd
import os
import yaml
import sys
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def inspect_read_types(sample, samplesheet):
    read1_from_sheet = samplesheet.at[sample, "read1"]
    read2_from_sheet = samplesheet.at[sample, "read2"]

    type = None

    if not pd.isna(read1_from_sheet) and Path(read1_from_sheet).exists():
        type = "SR"
        if not pd.isna(read2_from_sheet) and Path(read2_from_sheet).exists():
            type = "PR"
    elif not pd.isna(read2_from_sheet) and Path(read2_from_sheet).exists():
        logger.error("Read 1 does not exist but read 2 does, is your samplesheet corrupt?")
        sys.exit(0)

    return type


def inspect_samplesheet_assembly_path(sample, samplesheet):
    """
    Inspects the assembly path for a given sample from the samplesheet.

    Args:
        sample (str): Sample name.
        samplesheet (pd.DataFrame): DataFrame containing the samplesheet data.

    Returns:
        dict: Dictionary with sample as key and Path, None, or False as value.
    """
    # Reading assembly entry from samplesheet
    assembly_from_sheet = samplesheet.at[sample, "assembly"]
    
    # Handle if assembly is determined as NA
    assembly_path = Path(str(assembly_from_sheet))
    
    if pd.isna(assembly_from_sheet):
        logger.trace(f"No assembly provided for {sample}")
        path = {sample: None}
    elif assembly_path.exists(follow_symlinks = True):
        path = {sample: assembly_path}
    else:
        logger.warning(f"Failed to find {assembly_from_sheet}")
        path = {sample: False}

    return path


def determine_sample_configs(samplesheet, config_dir, ignore_assemblies):
    """
    Determines the sample configurations by reading config files and handling assemblies.

    Args:
        samplesheet (pd.DataFrame): DataFrame containing the samplesheet data.
        config_dir (str): Directory path containing configuration files.
        ignore_assemblies (bool): Whether to ignore assemblies or not.

    Returns:
        dict: Dictionary with sample names as keys and their configurations as values.
    """
    
    # Create a dict for sample names and dict files
    sample_configs = {}
    assemblers_unknown = []

    # Iterate samplesheet and pair samples with configurations
    for sample, cfg in zip(samplesheet.index, samplesheet["config"]):

        # Determine assembly paths
        assembly_path = inspect_samplesheet_assembly_path(sample, samplesheet)

        # Add sample assembly to list of unknown assembler, if assembly exists
        if isinstance(assembly_path.get(sample), Path) and not ignore_assemblies:
            assemblers_unknown.append(sample)

        # Deduce configuration file from samplesheet
        cfg_path = f"{config_dir}/{cfg}"

        # Handle missing configuration file
        if not os.path.isfile(cfg_path):
            available = sorted(p.name for p in Path(config_dir).glob("*.yaml"))
            logger.warning(
                f"Sample '{sample}' references config '{cfg}' which was not "
                f"found in {config_dir}. Available configs: "
                f"{', '.join(available) if available else '(none)'}."
            )
            cfg_path = None

            default_path = f"{config_dir}/default.yaml"

            # Ensure that default file exists and use it
            if os.path.exists(default_path):
                logger.warning(f"Falling back to default.yaml for sample '{sample}'.")
                cfg_path = default_path
            else:
                logger.error(
                    "Default configuration file is missing, "
                    "please recreate it to enable default analysis: "
                    f"{default_path}"
                )


        # Read sample configurations
        if cfg_path is not None:
            # print(f"Configuration file {cfg_path} found for {sample}")  # As log_debug
            with open(cfg_path, "r") as config_file:
                sample_cfg = yaml.safe_load(config_file)
                sample_configs[sample] = sample_cfg

            for mod, opts in sample_cfg.items():

                if not isinstance(opts, dict):
                    sample_configs[sample][mod] = {"enabled": opts}

                # Record read type
                sample_configs[sample][mod]["read_type"] = inspect_read_types(sample, samplesheet)

        # Warn user of no configuration is included
        else:
            print(f"No configuration file was specified for sample {sample}.\nSkipping!")
    
    # Ensure that there are indeed sample configurations
    if len(sample_configs) == 0:
        sys.exit("No sample configuration files found. Ensure that the `config` column of the samplesheet is correctly filled.")

    # Changing assembler values for samples with unknown assemblers
    for sample in assemblers_unknown:
        sample_cfg = sample_configs.get(sample)

        for mod, opts in sample_cfg.items():
            if not isinstance(opts, dict):
                continue
            elif "assembler" in opts.keys():
                sample_configs[sample][mod]["assembler"] = ["unkwnAssmbly"]

    return sample_configs
