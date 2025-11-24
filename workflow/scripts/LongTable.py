#!/usr/bin/env python3
import argparse
import os
import yaml
import pandas as pd
import re
from typing import List, Optional, Any, Dict, Tuple
from pandas.errors import EmptyDataError

# ---------- helpers ----------
def analysis_to_run_check(entry: Any) -> bool:
    """ check if a analysis_to_run is true to run that particular tool """
    return entry is True or isinstance(entry, dict)

def parse_option_species_value(options: str, keys: List[str]) -> Optional[str]:
    """
    parse tool specific config entries in options related to species or organism 
    as some tools use it for the output name
    """
    if not options:
        return None
    for k in keys:
        match = re.search(rf"{re.escape(k)}\s+'([^']+)'", options) # --species 'Escherichia coli'
        if match: return match.group(1)
        match = re.search(rf'{re.escape(k)}\s+"([^"]+)"', options) # --species "Escherichia coli"
        if match: return match.group(1)
        match = re.search(rf"{re.escape(k)}\s+([^\s]+)", options) # --species Escherichia coli
        if match: return match.group(1)
    return None

def space_to_underscore(text: Optional[str]) -> Optional[str]:
    """
    replace space in organism name with underscore to match output filename
    """
    return re.sub(r"\s+", "_", text.strip()) if text else None

def infer_read1(illumina_files: Optional[str]) -> Optional[str]:
    """
    some tools use the full name of the first read for names of the output files
    with the values extracted from the sample sheet
    """
    if not illumina_files:
        return None
    first = illumina_files.split(",")[0].strip()
    return os.path.basename(first)

def convert_to_list(value: Any) -> List[Optional[str]]:
    """Convert into a list."""
    if value is None:
        return [None]
    if isinstance(value, list):
        return value if value else [None]
    return [str(value)]

def value_is_na(x: Optional[str]) -> bool:
    """Return True if x is empty or an NA like value (NA..NA, nan, NA) as case-insensitive to remove it."""
    if x is None:
        return True
    s = str(x).strip().lower()
    return s in {"", "na", "n/a", "nan", "na..na"}

# ---------- file reading ----------
def Check_emptylines_header(path: str) -> Tuple[int, Optional[str]]:
    """ returns the number of non-empty lines, and the first non-empty line itself -
    which we use to check if a file has any real data (non-empty lines), 
    and whether it’s header-only (just one line of column names).
    """

    lines = 0
    first = None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            s = raw.strip("\n\r")
            if s.strip() == "":
                continue
            if first is None:
                first = s
            lines += 1
    return lines, first

def read_table_mlst(path: str) -> Optional[pd.DataFrame]:
    """
    FOR MLST ONLY: read a headerless single- or multi-line table.
    Accepts TSV or whitespace-delimited. Assigns generic col names.
    Skips zero-byte files; otherwise tries hard not to drop data.
    """
    try:
        if os.path.getsize(path) == 0:
            print(f"\t  Skipping empty file: {path}")
            return None
    except OSError:
        pass

    # Try TSV header=None
    try:
        df = pd.read_csv(path, sep="\t", header=None, dtype=str, engine="python")
        if df.shape[1] >= 1 and df.shape[0] >= 1:
            df.columns = [f"col{i+1}" for i in range(df.shape[1])]
            return df
    except Exception:
        pass

    # Try whitespace header=None
    try:
        df2 = pd.read_csv(path, delim_whitespace=True, header=None, dtype=str, engine="python")
        if df2.shape[1] >= 1 and df2.shape[0] >= 1:
            df2.columns = [f"col{i+1}" for i in range(df2.shape[1])]
            return df2
    except Exception:
        pass

    # Last resort: if there is exactly one non-empty line, split it
    try:
        nonempty_count, first = Check_emptylines_header(path)
        if nonempty_count == 1 and first is not None:
            row = first.split("\t") if ("\t" in first) else first.split()
            if len(row) >= 1:
                df1 = pd.DataFrame([row], columns=[f"col{i+1}" for i in range(len(row))])
                return df1
    except Exception:
        pass

    print(f"\t  Skipping unreadable headerless file: {path}")
    return None

def read_table_information(path: str) -> Optional[pd.DataFrame]:
    """
    General reader for non-MLST tools:
      1) TSV with header
      2) TSV header=None (to rescue headerless data)
      3) Whitespace-delimited header=None
      4) CSV with header
    Skips zero-byte and truly column-less files. Avoids including files
    that are header-only when parsed as TSV with header (0 rows).
    """
    # Skip zero-byte files
    try:
        if os.path.getsize(path) == 0:
            print(f"\t  Skipping empty file: {path}")
            return None
    except OSError:
        pass

    # Attempt 1: TSV with header
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        if df.shape[1] == 0:
            print(f"\t  Skipping file with no columns: {path}")
            return None
        if df.shape[0] > 0:
            return df
        # If zero rows, check if header-only
        nonempty_count, _ = Check_emptylines_header(path)
        if nonempty_count <= 1:
            print(f"\t  Skipping header-only file: {path}")
            return None
        # else: fall through to header=None
    except EmptyDataError:
        # Try fallbacks
        pass
    except Exception:
        # Try fallbacks
        pass

    # Attempt 2: TSV header=None (allows single-line data)
    try:
        df2 = pd.read_csv(path, sep="\t", header=None, dtype=str, engine="python")
        if df2.shape[1] >= 1 and df2.shape[0] >= 1:
            df2.columns = [f"col{i+1}" for i in range(df2.shape[1])]
            return df2
    except Exception:
        pass

    # Attempt 3: whitespace-delimited header=None
    try:
        df3 = pd.read_csv(path, delim_whitespace=True, header=None, dtype=str, engine="python")
        if df3.shape[1] >= 1 and df3.shape[0] >= 1:
            df3.columns = [f"col{i+1}" for i in range(df3.shape[1])]
            return df3
    except Exception:
        pass

    # Attempt 4: CSV with header
    try:
        df4 = pd.read_csv(path, sep=",", dtype=str, low_memory=False)
        if df4.shape[1] >= 1 and df4.shape[0] >= 1:
            return df4
    except Exception:
        pass

    print(f"\t  Skipping unreadable file after fallbacks: {path}")
    return None

# ---------- pattern resolution ----------
def resolve_patterns(patterns: List[str], tool_cfg: dict, sample_row: pd.Series) -> List[str]:
    asm_vals = convert_to_list(tool_cfg["assemblers"]) if "assemblers" in tool_cfg else [None]
    db_vals  = convert_to_list(tool_cfg["database"])   if "database"   in tool_cfg else [None]

    if "organism" in tool_cfg:
        organism_text = str(tool_cfg["organism"])
    else:
        options = tool_cfg["options"] if "options" in tool_cfg else ""
        organism_text = parse_option_species_value(options, ["--species", "--organism"]) or None
    organism_us = space_to_underscore(organism_text)

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
    seen = set()
    uniq = []
    for r in resolved:
        if r not in seen:
            seen.add(r)
            uniq.append(r)
    return uniq

# ---------- long_table_creation ----------
def long_table_creation(df: pd.DataFrame,
                            tool: str,
                            sample_id_fallback: str,
                            organism_fallback: Optional[str]) -> pd.DataFrame:
    df = df.copy().astype(str)

    # Case-insensitive lookup of existing id columns
    colmap = {c.lower(): c for c in df.columns}
    sid_col = colmap.get("sample_name")
    org_col = colmap.get("organism")

    # Add missing id columns in one shot
    to_add = {}
    if not sid_col:
        to_add["sample_name"] = sample_id_fallback
    if not org_col:
        to_add["organism"] = organism_fallback if organism_fallback is not None else ""
    if to_add:
        df = df.assign(**to_add)

    # Recompute names in case we added columns
    colmap = {c.lower(): c for c in df.columns}
    sid_col = colmap.get("sample_name", "sample_name")
    org_col = colmap.get("organism", "organism")

    id_vars = [sid_col, org_col]
    value_vars = [c for c in df.columns if c not in id_vars]

    long = df.melt(id_vars=id_vars, value_vars=value_vars,
                   var_name="field", value_name="value")

    long = (
        long.rename(columns={sid_col: "sample_name", org_col: "organism"})
            .assign(tool=tool, filename="")
            [["sample_name", "tool", "filename", "organism", "field", "value"]]
    )
    return long

# ---------- core work moved here ----------
def save_long_table(
    samplesheet_path: str,
    catalogue_path: str,
    config_species_root: str,
    output_folder: str,
    na_filter: bool,
    extend: bool,
) -> None:

    if not os.path.exists(samplesheet_path):
        raise FileNotFoundError(f"Samplesheet not found: {samplesheet_path}")
    if not os.path.exists(catalogue_path):
        raise FileNotFoundError(f"Catalogue not found: {catalogue_path}")

    os.makedirs(output_folder, exist_ok=True)

    samplesheet = pd.read_csv(samplesheet_path, sep="\t")
    print(f"Using samplesheet: {samplesheet_path}")
    print(f"\nSuccessfully read samplesheet with {samplesheet.shape[0]} rows and {samplesheet.shape[1]} columns")
    print(samplesheet.head(), "\n")

    with open(catalogue_path, "r") as f:
        catalogue = yaml.safe_load(f)
    print(f"Using catalogue: {catalogue_path}")
    print(f"\nSuccessfully read results catalogue with {len(catalogue)} entries.\n")

    combined_frames: List[pd.DataFrame] = []  # used only when extend=True
    written: List[str] = []

    for _, row in samplesheet.iterrows():
        sample = str(row["sample_name"])
        species_cfg_file = str(row["config"])
        species_cfg_path = os.path.join(config_species_root, species_cfg_file)
        if not os.path.exists(species_cfg_path):
            raise FileNotFoundError(f"Species config not found for sample {sample}: {species_cfg_path}")

        with open(species_cfg_path, "r") as f:
            species_cfg = yaml.safe_load(f)

        analyses = species_cfg #load in all entries, before it was analysis_to_run
        enabled_tools = [tool for tool, entry in analyses.items() if analysis_to_run_check(entry)]
        relevant_tools = [tool for tool in enabled_tools if tool in catalogue]

        print(f"=== Sample: {sample} | species config: {species_cfg_file} ===")
        if not relevant_tools:
            print("(No relevant catalogue entries)\n")
        else:
            print("Relevant catalogue entries:")

        # per-tool organism fallback (pretty)
        tool_org_fallbacks: Dict[str, Optional[str]] = {}
        for tool in relevant_tools:
            entry = analyses[tool]
            if isinstance(entry, dict) and "organism" in entry:
                tool_org_fallbacks[tool] = str(entry["organism"]).replace("_", " ")
            else:
                options = entry["options"] if (isinstance(entry, dict) and "options" in entry) else ""
                tool_org_fallbacks[tool] = parse_option_species_value(options, ["--species", "--organism"])

        sample_frames: List[pd.DataFrame] = []
        for tool in relevant_tools:
            tool_cfg = analyses[tool] if isinstance(analyses[tool], dict) else {}
            cat_val = catalogue[tool]
            patterns = cat_val if isinstance(cat_val, list) else [cat_val]
            rels = resolve_patterns(patterns, tool_cfg, row)
            base = os.path.join(output_folder, sample, tool)
            files = [os.path.join(base, r) for r in rels]

            # Print resolved paths
            if len(files) == 1:
                print(f"\t-{tool}: {files[0]}")
            else:
                print(f"\t-{tool}:")
                for fpath in files:
                    print(f"\t\t• {fpath}")

            # Read files
            for fp in files:
                if not os.path.isfile(fp):
                    continue

                if tool == "mlst":
                    df = read_table_mlst(fp)
                else:
                    df = read_table_information(fp)
                if df is None:
                    continue

                longtable = long_table_creation(
                    df=df,
                    tool=tool,
                    sample_id_fallback=sample,
                    organism_fallback=tool_org_fallbacks[tool]
                )
                longtable["filename"] = os.path.basename(fp)
                sample_frames.append(longtable)

        print("")

        # Build per-sample table
        if sample_frames:
            out_df = pd.concat(sample_frames, ignore_index=True)
        else:
            out_df = pd.DataFrame(columns=["sample_name", "tool", "filename", "organism", "field", "value"])

        if na_filter and not out_df.empty:
            before = len(out_df)
            out_df = out_df[~out_df["value"].apply(value_is_na)].copy()
            removed = before - len(out_df)
            if removed > 0:
                print(f"  (Filtered {removed} NA-like rows for {sample} due to --na)")

        if extend:
            # collect for one combined write/append later
            combined_frames.append(out_df)
        else:
            # write per-sample file
            out_path = os.path.join(output_folder, f"{sample}/{sample}_longtable.tsv")
            out_df.to_csv(out_path, sep="\t", index=False)
            written.append(out_path)

    if extend:
        # Write/append the single combined file
        combined = pd.concat(combined_frames, ignore_index=True) if combined_frames else \
                   pd.DataFrame(columns=["sample_name", "tool", "filename", "organism", "field", "value"])
        out_path = os.path.join(output_folder, "all_results.tsv")
        exists_before = os.path.exists(out_path)

        if exists_before:
            combined.to_csv(out_path, sep="\t", index=False, mode="a", header=False)
            print(f"Appended {len(combined)} rows to existing file: {out_path}")
        else:
            combined.to_csv(out_path, sep="\t", index=False)
            print(f"Wrote {len(combined)} rows to: {out_path}")
        written.append(out_path)

    print("Long-table files written:")
    for path in written:
        try:
            n = sum(1 for _ in open(path)) - 1
            print(f"  - {path}  ({max(n,0)} rows)")
        except Exception:
            print(f"  - {path}")


# ---------- main ----------
def main():
    parser = argparse.ArgumentParser(
        description="Resolve outputs per sample and emit per-sample or combined long-table TSVs."
    )
    parser.add_argument("--config", required=True, help="Path to main config YAML")
    parser.add_argument("--sample_sheet", help="Optional path to samplesheet (TSV)")
    parser.add_argument("--catalogue", help="Optional path to results catalogue YAML")
    parser.add_argument(
        "-n", "--na",
        action="store_true",
        help="If set, drop rows where the 'value' field is empty or NA-like (NA, N/A, NaN, NA..NA)."
    )
    parser.add_argument(
        "-e", "--extend",
        action="store_true",
        help="If set, append all samples into {output_folder}/Longtable.txt (create if missing). Without this flag, write per-sample files."
    )
    args = parser.parse_args()

    # Read config
    if not os.path.exists(args.config):
        raise FileNotFoundError(f"Config file not found: {args.config}")
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)
    print(f"Successfully read config file: {args.config}\n")

    samplesheet_path = args.sample_sheet or config["samplesheet"]
    catalogue_path = args.catalogue or config["result_catalogue"]
    config_species_root = config["config_species"]
    output_folder = config["output_folder"]

    save_long_table(
        samplesheet_path=samplesheet_path,
        catalogue_path=catalogue_path,
        config_species_root=config_species_root,
        output_folder=output_folder,
        na_filter=args.na,
        extend=args.extend,
    )

if __name__ == "__main__":
    main()