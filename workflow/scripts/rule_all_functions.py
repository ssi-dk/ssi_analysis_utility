import pandas as pd
import os
import yaml

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
