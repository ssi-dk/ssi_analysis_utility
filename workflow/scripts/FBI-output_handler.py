import pandas as pd
import argparse
import os


# Settings
fliflag = "false"

# Thresholds for various genes
thresholds = {
    "stx": [98, 98],
    "wzx": [98, 98],
    "wzy": [98, 98],
    "wzt": [98, 98],
    "wzm": [98, 98],
    "fliC": [90, 90],
    "fli": [90, 90],
    "eae": [95, 95],
    "ehxA": [95, 95],
    "other": [98, 98]
}

stxbadthreshold = [30, 90]


def get_threshold(template_name, thresholds):
    """Returns the appropriate thresholds based on the template name."""
    for key in thresholds:
        if key in template_name:
            return thresholds[key]
    return thresholds["other"]


def process_res_file(res_file_path):
    """Process the .res file and filter based on thresholds."""
    try:
        df = pd.read_csv(res_file_path, sep="\t")
    except FileNotFoundError:
        print(f"Error: File {res_file_path} not found.")
        return None
    except pd.errors.EmptyDataError:
        print(f"Error: File {res_file_path} is empty or not properly formatted.")
        return None

    # Apply filtering based on threshold for Template_Coverage and Query_Identity
    df["threshold"] = df["#Template"].apply(lambda x: get_threshold(x, thresholds))
    df_filtered = df[
        (df["Template_Coverage"] >= df["threshold"].apply(lambda x: x[0])) &
        (df["Query_Identity"] >= df["threshold"].apply(lambda x: x[1]))
    ]
    
    return df_filtered

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='FBI-output_handler.py')
    parser.add_argument("-i", "--inputfile", 
                        help='Path to .res file', 
                        required=True)
    parser.add_argument("-s", "--sampleID", 
                        help='Sample Name', 
                        required=True)
    parser.add_argument("-st","--ST_list" ,
                        help="ST list", 
                        default=['ST:NA', 'NA:NA', 'NA:NA', 'NA:NA', 'NA:NA', 'NA:NA', 'NA:NA', 'NA:NA']) # ['ST:11', 'adk:12', 'fumC:12', 'gyrB:8', 'icd:12', 'mdh:15', 'purA:2', 'recA:2']
    args = parser.parse_args()

    # Define file path
    res_file_path = os.path.abspath(args.inputfile)

    # Process the .res file and get the filtered DataFrame
    filtered_df = process_res_file(res_file_path)    



    d = {}
    d["SampleID"] = args.sampleID
    
    # Create a dictionary to map gene names to output keys
    gene_map = {
        "wzx": "wzx",
        "wzy": "wzy",
        "wzt": "wzt",
        "wzm": "wzm",
        "fliC": "fliC",
        "fli": "fli",
        "eae": "eae",
        "ehxA": "ehxA"
    }

    # Iterate over the filtered DataFrame
    for template in filtered_df["#Template"]:
        gene_list = template.split("__")
        gene_name = gene_list[1]

        # Check if the gene name matches any in the map
        if gene_name in gene_map:
            # If the value is 'Positive', update accordingly
            if gene_name in ["eae", "ehxA"]:
                d[gene_map[gene_name]] = ["Positive"]
            else:
                d[gene_map[gene_name]] = [gene_list[2]]
        elif gene_name.startswith("stx"):
            d["stx"] = [gene_list[2]]
        elif gene_name not in thresholds:
            d["Other"] = [gene_list[2]]
        

    d["ST"] = args.ST_list.pop(0).split(":")[1]
    d["STstring"] = ",".join(args.ST_list)

    df1 = pd.DataFrame.from_dict(d)
    print(df1)