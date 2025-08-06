import pandas as pd
import os
import yaml
import warnings

# Tools that produce .sam and are used downstream in samtools, bcftools, etc.
Tools_with_sam_output = ["kmeraligner"]

# Assemblers that can be used for tools like AMRFinder or Repeat_identifier
Assembly_tools = ["spades", "skesa"]

def load_inputs(config):
    """
    Load inputs and data from the config and CSV samplesheet.

    Returns:
        CSV_DATA (pd.DataFrame): Raw parsed sample sheet
        SAMPLES (list): List of sample names
        sample_to_illumina (dict): Map sample -> [read1, read2]
        sample_to_nanopore (dict): Map sample -> nanopore file
        sample_to_assembly_file (dict): Map sample -> assembly FASTA
        sample_to_organism (dict): Map sample -> normalized species name
    """

    # read the csv file
    csv_path = config['input_manager']['path']
    csv_data = pd.read_csv(csv_path, sep='\t')

    # define relevant paths 
    output_folder = config['input_manager']['out_folder']
    database_path = config['input_manager']['database_path']

    samples = csv_data['sample_name'].tolist()

    # Create a dictionary for mapping sampleID to file paths
    sample_to_illumina = {
        row['sample_name']: row['Illumina_read_files'].split(',')
        for _, row in csv_data.iterrows()
    }

    sample_to_nanopore = {
        row['sample_name']: row['Nanopore_read_file']
        for _, row in csv_data.iterrows()
    }

    sample_to_assembly_file = {
        row['sample_name']: row['assembly_file']
        for _, row in csv_data.iterrows()
    }

    # Load species_name_map from config file
    species_name_map = config.get("species_map", {})
    
    # Map each sample species (from sample sheet information) to the normalized species name used for the Species specific config
    sample_to_organism = {
        row['sample_name']: species_name_map.get(row['organism'], row['organism'])
        for _, row in csv_data.iterrows()
    }

    # Debug check if species have not properly been mapped
    unmapped_species = {
        row['organism'] for _, row in csv_data.iterrows()
        if row['organism'] not in species_name_map
    }
    if unmapped_species:
        warnings.warn(
            f"The following species are not mapped in config['species_map']: {unmapped_species}",
            UserWarning
        )

    return (
        csv_data,
        samples,
        sample_to_illumina,
        sample_to_nanopore,
        sample_to_assembly_file,
        sample_to_organism,
        output_folder,
        database_path
    )

def load_species_configs(species_config_path, sample_to_organism):
    """
    Load per-species config YAML files.

    Args:
        species_config_path (str): Folder containing species .yaml files
        sample_to_organism (dict): Map sample -> normalized species name

    Returns:
        dict: Map normalized species name -> config dictionary
    """
    species_configs = {}

    for species in set(sample_to_organism.values()):
        species_file = os.path.join(species_config_path, f"{species}.yaml")
        if os.path.exists(species_file):
            with open(species_file, "r") as f:
                species_configs[species] = yaml.load(f, Loader=yaml.Loader)
        else:
            warnings.warn(
                f"Warning: Configuration file {species_file} not found. Skipping analyses for {species}.",
                UserWarning
            )
            species_configs[species] = {}  # fallback: empty config

    return species_configs

def is_tool_enabled(sample, tool_name, sample_to_organism, species_configs):
    """
    Determine whether a given tool is enabled for a sample
    based on species-specific YAML config.

    Args:
        sample (str): Sample name
        tool_name (str): Tool name from config['analyses_to_run']
        sample_to_organism (dict)
        species_configs (dict)

    Returns:
        bool: True if enabled
    """
    organism = sample_to_organism.get(sample)
    species_config = species_configs.get(organism, {})
    tool_config = species_config.get("analyses_to_run", {}).get(tool_name)

    return tool_config is not None and tool_config.get("status", False) is True
