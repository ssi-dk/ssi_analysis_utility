#!/usr/bin/env python3
import argparse
import os
import yaml
import pandas as pd

def is_enabled(entry):
    """Return True if an analysis entry is enabled."""
    return entry is True or isinstance(entry, dict)

def main():
    parser = argparse.ArgumentParser(
        description="Read config, samplesheet, and catalogue; print relevant catalogue entries per sample."
    )
    parser.add_argument("--config", required=True, help="Path to main config YAML")
    parser.add_argument("--sample_sheet", help="Optional path to samplesheet (TSV)")
    parser.add_argument("--catalogue", help="Optional path to results catalogue YAML")
    args = parser.parse_args()

    # --- Check that config file exists ---
    if not os.path.exists(args.config):
        raise FileNotFoundError(f"Config file not found: {args.config}")

    # --- Read config ---
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)
    print(f"Successfully read config file: {args.config}\n")

    # --- Print Input Manager section ---
    input_manager = config["input_manager"]
    print("Input manager section:")
    for key, value in input_manager.items():
        print(f"  {key}: {value}")
    print()

    # --- Resolve inputs ---
    samplesheet_path = args.sample_sheet or input_manager["samplesheet"]
    catalogue_path   = args.catalogue or input_manager["result_catalogue"]
    config_species_root = input_manager["config_species"]

    # --- Check files exist ---
    if not samplesheet_path or not os.path.exists(samplesheet_path):
        raise FileNotFoundError(f" Samplesheet not found: {samplesheet_path}")
    if not catalogue_path or not os.path.exists(catalogue_path):
        raise FileNotFoundError(f" Catalogue not found: {catalogue_path}")

    # --- Read samplesheet ---
    samplesheet = pd.read_csv(samplesheet_path, sep="\t")
    print(f"Using samplesheet: {samplesheet_path}")
    print(f"\n Successfully read samplesheet with {samplesheet.shape[0]} rows and {samplesheet.shape[1]} columns")
    print(samplesheet.head(), "\n")

    # --- Read catalogue ---
    with open(catalogue_path, "r") as f:
        catalogue = yaml.safe_load(f)
    print(f"Using catalogue: {catalogue_path}")
    print(f"\n Successfully read results catalogue with {len(catalogue)} entries.\n")

    # --- For each sample: load species config and print relevant catalogue entries ---
    for _, row in samplesheet.iterrows():
        sample = str(row["sample_name"])
        species_cfg_file = str(row["config"])
        species_cfg_path = os.path.join(config_species_root, species_cfg_file)

        if not os.path.exists(species_cfg_path):
            raise FileNotFoundError(f"Species config not found for sample {sample}: {species_cfg_path}")

        with open(species_cfg_path, "r") as f:
            species_cfg = yaml.safe_load(f) or {}

        analyses = species_cfg["analyses_to_run"]
        enabled_tools = [tool for tool, entry in analyses.items() if is_enabled(entry)]
        relevant = [tool for tool in enabled_tools if tool in catalogue]

        print(f"=== Sample: {sample} | species config: {species_cfg_file} ===")
        if not relevant:
            print("  (No relevant catalogue entries for enabled analyses)\n")
            continue

        print("  Relevant catalogue entries:")
        for tool in relevant:
            cat_val = catalogue[tool]
            if isinstance(cat_val, list):
                print(f"  - {tool}:")
                for item in cat_val:
                    print(f"      • {item}")
            else:
                print(f"  - {tool}: {cat_val}")
        print("")
if __name__ == "__main__":
    main()
