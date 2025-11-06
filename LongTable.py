#!/usr/bin/env python3
import argparse
import os
import yaml
import pandas as pd
import re
from typing import List, Optional, Any

def is_enabled(entry: Any) -> bool:
    """Enabled if True or a dict with settings/options."""
    return entry is True or isinstance(entry, dict)

def parse_option_value(options: str, keys: List[str]) -> Optional[str]:
    """Extract the value for any of keys (e.g., --species/--organism) from an options string."""
    if not options:
        return None
    for k in keys:
        # support 'single quotes', "double quotes", and unquoted single token
        m = re.search(rf"{re.escape(k)}\s+'([^']+)'", options)
        if m: return m.group(1)
        m = re.search(rf'{re.escape(k)}\s+"([^"]+)"', options)
        if m: return m.group(1)
        m = re.search(rf"{re.escape(k)}\s+([^\s]+)", options)
        if m: return m.group(1)
    return None

def to_underscore(text: Optional[str]) -> Optional[str]:
    # resfinder defines in its pointfinder the species in the output filename
    # examples/Results/SRR26205262/resfinder/pointfinder_kma/kma_Salmonella_enterica_SRR26205262_1.fastq.gz.res
    return re.sub(r"\s+", "_", text.strip()) if text else None 

def infer_read1(illumina_files: Optional[str]) -> Optional[str]:
    """Return the full first read filename (including extensions)."""
    # resfinder defines in its pointfinder the full read1 file name with extension in the output filename
    # examples/Results/SRR26205262/resfinder/pointfinder_kma/kma_Salmonella_enterica_SRR26205262_1.fastq.gz.res
    if not illumina_files:
        return None
    first = illumina_files.split(",")[0].strip()
    return os.path.basename(first)

def listify(value: Any) -> List[Optional[str]]:
    """Normalize a scalar or list into a list (None -> [None])."""
    if value is None:
        return [None]
    if isinstance(value, list):
        return value if value else [None]
    return [str(value)]

def resolve_patterns(patterns: List[str], tool_cfg: dict, sample_row: pd.Series) -> List[str]:
    """
    Substitute placeholders in catalogue patterns using per-tool config + sample row.
    Handles {assemblers}, {database}, {organism}, {read1}. Fans out all combinations.
    """
    asm_vals = listify(tool_cfg["assemblers"]) if "assemblers" in tool_cfg else [None]
    db_vals  = listify(tool_cfg["database"])   if "database"   in tool_cfg else [None]

    # organism: prefer explicit key, else derive from options
    if "organism" in tool_cfg:
        organism = str(tool_cfg["organism"])
    else:
        options = tool_cfg["options"] if "options" in tool_cfg else ""
        organism = parse_option_value(options, ["--species", "--organism"]) or None
    organism_us = to_underscore(organism)

    read1 = infer_read1(sample_row["Illumina_read_files"]) if "Illumina_read_files" in sample_row else None

    resolved: List[str] = []
    for pat in patterns:
        for asm in asm_vals:
            for db in db_vals:
                out = pat
                if "{assemblers}" in out and asm is not None:
                    out = out.replace("{assemblers}", str(asm))
                if "{database}" in out and db is not None:
                    out = out.replace("{database}", str(db))
                if "{organism}" in out and organism_us is not None:
                    out = out.replace("{organism}", organism_us)
                if "{read1}" in out and read1 is not None:
                    out = out.replace("{read1}", read1)
                resolved.append(out)
    # Deduplicate while preserving order
    seen = set()
    uniq = []
    for r in resolved:
        if r not in seen:
            seen.add(r)
            uniq.append(r)
    return uniq

def main():
    parser = argparse.ArgumentParser(
        description="Read config, samplesheet, and catalogue; print resolved catalogue entries per sample (with full output paths)."
    )
    parser.add_argument("--config", required=True, help="Path to main config YAML")
    parser.add_argument("--sample_sheet", help="Optional path to samplesheet (TSV)")
    parser.add_argument("--catalogue", help="Optional path to results catalogue YAML")
    args = parser.parse_args()

    # --- Check & read main config (hard brackets everywhere) ---
    if not os.path.exists(args.config):
        raise FileNotFoundError(f"Config file not found: {args.config}")
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)
    print(f"Successfully read config file: {args.config}\n")

    input_manager = config["input_manager"]
    print("Input manager section:")
    for key, value in input_manager.items():
        print(f"  {key}: {value}")
    print()

    # Resolve paths
    samplesheet_path     = args.sample_sheet if args.sample_sheet else input_manager["samplesheet"]
    catalogue_path       = args.catalogue    if args.catalogue    else input_manager["result_catalogue"]
    config_species_root  = input_manager["config_species"]
    output_folder        = input_manager["output_folder"]

    # Existence checks
    if not os.path.exists(samplesheet_path):
        raise FileNotFoundError(f"Samplesheet not found: {samplesheet_path}")
    if not os.path.exists(catalogue_path):
        raise FileNotFoundError(f"Catalogue not found: {catalogue_path}")

    # Load data
    samplesheet = pd.read_csv(samplesheet_path, sep="\t")
    print(f"Using samplesheet: {samplesheet_path}")
    print(f"\nSuccessfully read samplesheet with {samplesheet.shape[0]} rows and {samplesheet.shape[1]} columns")
    print(samplesheet.head(), "\n")

    with open(catalogue_path, "r") as f:
        catalogue = yaml.safe_load(f)
    print(f"Using catalogue: {catalogue_path}")
    print(f"\nSuccessfully read results catalogue with {len(catalogue)} entries.\n")

    # Per-sample processing (now printing full paths under output_folder/sample/tool/…)
    for _, row in samplesheet.iterrows():
        sample = str(row["sample_name"])
        species_cfg_file = str(row["config"])
        species_cfg_path = os.path.join(config_species_root, species_cfg_file)
        if not os.path.exists(species_cfg_path):
            raise FileNotFoundError(f"Species config not found for sample {sample}: {species_cfg_path}")

        with open(species_cfg_path, "r") as f:
            species_cfg = yaml.safe_load(f)

        analyses = species_cfg["analyses_to_run"]
        enabled_tools = [tool for tool, entry in analyses.items() if is_enabled(entry)]
        relevant_tools = [tool for tool in enabled_tools if tool in catalogue]

        print(f"=== Sample: {sample} | species config: {species_cfg_file} ===")
        if not relevant_tools:
            print("  (No relevant catalogue entries for enabled analyses)\n")
            continue

        print("  Relevant catalogue entries:")
        for tool in relevant_tools:
            tool_cfg = analyses[tool] if isinstance(analyses[tool], dict) else {}
            cat_val = catalogue[tool]
            patterns = cat_val if isinstance(cat_val, list) else [cat_val]
            resolved = resolve_patterns(patterns, tool_cfg, row)

            # Join with output_folder/sample/tool
            base = os.path.join(output_folder, sample, tool)
            full_paths = [os.path.join(base, r) for r in resolved]

            if len(full_paths) == 1:
                print(f"\t-{tool}: {full_paths[0]}")
            else:
                print(f"\t-{tool}:")
                for p in full_paths:
                    print(f"\t\t• {p}")
        print("")

if __name__ == "__main__":
    main()
