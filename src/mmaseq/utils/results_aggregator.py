#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from collections import defaultdict

def unpivot_results(sample, module, file, results):
    """
    Converts a wide-format results DataFrame to a long-format DataFrame.

    Args:
        sample (str): Sample name.
        module (str): Module name.
        file (Path): Path to the results file.
        results (pd.DataFrame): Wide-format results DataFrame.

    Returns:
        pd.DataFrame: Long-format DataFrame with columns Sample, Module, File, Row, Column, Value.
    """
    # Convert index to row number column starting from row 1
    results.index += 1
    results = results.reset_index(names = "Row")

    # Generate long list format of results file
    results_long = results.melt(
        id_vars = "Row",
        var_name = "Column",
        value_name = "Value"
        )

    # add columns in a single assignment (faster than multiple insert calls)
    results_long[["Sample", "Module", "File"]] = [sample, module, file.name]
    
    # if you want a specific column order:
    cols = ["Sample", "Module", "File", "Row", "Column", "Value"]
    results_long = results_long[cols]

    return results_long



def generate_long_results(all_result_files):
    """
    Generates a concatenated long-format DataFrame from all result files.

    Args:
        all_result_files (defaultdict): Nested dictionary with sample -> module -> list of Path objects.

    Returns:
        pd.DataFrame: Concatenated long-format DataFrame from all result files.
    """
    all_sample_results = list()

    for sample, modules in all_result_files.items():

        for mod, files in modules.items():
            for file in files:
                try:
                    sample_results = pd.read_csv(file, sep = "\t", index_col = False)
                except pd.errors.EmptyDataError:
                    print(f"Results file {file.name} for {sample} is empty. Skipping!")
                    continue

                # Determine whether the long table format is already observed
                sample_long = sample_results
                required_columns = {
                    "Sample", "Module", "File", "Row", "Column", "Value"
                }
                if not required_columns.issubset(sample_results.columns):
                    sample_long = unpivot_results(sample, mod, file, sample_results)

                all_sample_results.append(sample_long)

    return pd.concat(all_sample_results, ignore_index = True)
