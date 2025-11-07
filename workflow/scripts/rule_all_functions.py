import pandas as pd
import os
import yaml
from typing import Dict, Tuple, Set, List
import warnings


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
    'examples/Results/ERR142064/deletion_identifier/Cdiff_Toxins.done', 
    'examples/Results/ERR142064/kma_filter/Cdiff_Toxins.done', 
    'examples/Results/ERR142064/mlst/spades.done', 
    'examples/Results/ERR142064/snp_identifier/Cdiff_Toxins.done', 
    'examples/Results/ERR3528110/amrfinder/shovill.done', 
    'examples/Results/ERR3528110/chtyper/fumCH_db.done', 
    'examples/Results/ERR3528110/custom_blaster/shovill_OXAndm.done',
    'examples/Results/ERR3528110/kma_filter/ecoligenes.done', 
    'examples/Results/ERR3528110/mlst/shovill.done', 
    'examples/Results/ERR3528110/plasmidfinder/plasmidfinder.done',
    'examples/Results/ERR3528110/resfinder/resfinder.done',
    'examples/Results/ERR3528110/serotypefinder/serotypefinder.done', 
    'examples/Results/ERR3528110/virulencefinder/virulencefinder.done', 
    'examples/Results/SRR25448586/amrfinder/shovill.done', 
    'examples/Results/SRR25448586/meningotype/shovill.done', 
    'examples/Results/SRR25448586/mlst/shovill.done', 
    'examples/Results/SRR25448586/plasmidfinder/plasmidfinder.done', 
    'examples/Results/SRR25448586/serotypefinder/serotypefinder.done', 
    'examples/Results/SRR25448586/virulencefinder/virulencefinder.done', 
    'examples/Results/SRR26205262/amrfinder/shovill.done', 
    'examples/Results/SRR26205262/mlst/shovill.done', 
    'examples/Results/SRR26205262/resfinder/resfinder.done', 
    'examples/Results/SRR26205262/seqsero2/seqsero2.done', 
    'examples/Results/SRR26205262/sistr/shovill.done', 
    'examples/Results/SRR4046826/amrfinder/shovill.done', 
    'examples/Results/SRR4046826/kleborate/shovill.done', 
    'examples/Results/SRR4046826/mlst/shovill.done', 
    'examples/Results/SRR4046826/plasmidfinder/plasmidfinder.done', 
    'examples/Results/SRR4046826/resfinder/resfinder.done', 
    'examples/Results/SRR4046826/serotypefinder/serotypefinder.done', 
    'examples/Results/SRR4046826/virulencefinder/virulencefinder.done']
    """
    return sorted(results)
