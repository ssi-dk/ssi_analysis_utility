import pandas as pd
import os
import yaml
from typing import Dict, Tuple, Set, List
import warnings


def resolve_env(envs_path,
                tool_name):
    """
    Return the conda environment path.
    If it doesn't exist, fall back to the workflow/envs/<tool_name>.yaml.
    """
    # Get path
    pipeline_dir = os.getcwd()
    print_path = False
    
    # Build path to env
    full_path = os.path.join(envs_path, tool_name)
    if os.path.exists(full_path):
        return tool_name
    else:
        # Revert to yaml if the env is not available
        yaml_path = "%s/workflow/envs/%s.yaml" % (pipeline_dir, tool_name)
        if os.path.exists(yaml_path):
            print(f"{tool_name}: Environment not found. Reverting to Snakemake installation. Disable this warning by setting envs_location to None.")
        else:
            # No yaml for the tool
            raise ValueError(f"No valid conda environment or YAML found for '{tool_name}'.")

    return yaml_path



def sample_read_map(
    samplesheet: pd.DataFrame,
    sample_col: str = "sample_name",
    read_col: str = "Illumina_read_files",
    nanopore_col: str = "Nanopore_read_file",
    assembly_col: str = "assembly_file"
) -> Tuple[Dict[str, List[str]], Dict[str, str], Dict[str, str]]:

    #{'ERR3528110': ['examples/Dataset/reads/ERR3528110_1.fastq.gz', 'examples/Dataset/reads/ERR3528110_2.fastq.gz']}
    sample_to_illumina = {
        row[sample_col]: row[read_col].split(',')
        for idx, row in samplesheet.iterrows()
    }
    
    sample_to_nanopore = {
        row[sample_col]: row[nanopore_col] 
        for idx, row in samplesheet.iterrows()
    }

    sample_to_assembly_file = {
        row[sample_col]: row[assembly_col] 
        for idx, row in samplesheet.iterrows()
    }

    return sample_to_illumina, sample_to_nanopore, sample_to_assembly_file


def sample_map(
    species_name_map: Dict[str, str],
    samplesheet: pd.DataFrame,
    sample_col: str = "sample_name",
    organism_col: str = "organism"
) -> Tuple[Dict[str, str], Set[str]]:
    """
    Map each sample to a normalized species name and return any unmapped raw species.
      - sample_to_organism: {sample: normalized_species}
      - unmapped_species: set of raw species not present in species_name_map
    """
    sample_to_organism = {
        row[sample_col]: species_name_map.get(row[organism_col], row[organism_col])
        for idx, row in samplesheet.iterrows()
    }

    unmapped_species = {
        row[organism_col] for idx, row in samplesheet.iterrows()
        if row[organism_col] not in species_name_map
    }

    return sample_to_organism, unmapped_species


def sample_to_species_config(
    sample_to_organism : dict,
    species_config_path : str = 'workflow/configs_species/'
    ) -> Dict:
    species_configs = {}
    for species in set(sample_to_organism.values()):
       
        species_config_file = os.path.join(species_config_path, f"{species}.yaml")

        if os.path.exists(species_config_file):
            with open(species_config_file, "r") as species_config_file:
                species_configs[species] = yaml.load(species_config_file,Loader=yaml.Loader)
        else: # means the config is not found / available
            warnings.warn("Warning: Configuration file %s not found. Skipping analyses for %s." % (species_config_file, species), UserWarning)
            species_configs[species] = {}  # Empty config to ensure skipping analyses

    return species_configs

def list_results(samplesheet, species_configs, output_folder, species_name_map):
    results = set()

    for _, row in samplesheet.iterrows():
        sample = row["sample_name"]
        organism_raw = row["organism"]
        organism_norm = species_name_map.get(organism_raw, organism_raw)

        config = species_configs.get(organism_norm, {}) or {}
        analyses = config.get("analyses_to_run", {}) or {}

        for analysis_name, settings in analyses.items():
            # e.g. PlasmidFinder: true - so where no db currently is needed
            if not isinstance(settings, dict):
                results.add(f"{output_folder}/{sample}/{analysis_name}/{analysis_name}.done")
                continue

            assemblers = settings.get("assemblers")
            databases = settings.get("database")

            # ensure assemblers are a list of strings
            if isinstance(assemblers, str):
                assemblers = [assemblers]

            if databases is None:
                 # no database case
                if assemblers:
                    for tool in assemblers:
                        results.add(f"{output_folder}/{sample}/{analysis_name}/{tool}.done")
                else:
                    results.add(f"{output_folder}/{sample}/{analysis_name}/{analysis_name}.done")
            else:
                # database case (can be single string or list of strings)
                if isinstance(databases, str):
                    databases = [databases]
                if assemblers:
                    for tool in assemblers:
                        for database in databases:
                            results.add(f"{output_folder}/{sample}/{analysis_name}/{tool}_{database}.done")
                else:
                    for database in databases:
                        results.add(f"{output_folder}/{sample}/{analysis_name}/{database}.done")

    return sorted(results)
