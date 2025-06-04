from typing import List, Dict

# Thresholds for E. coli [Template_Coverage, Query_Identity]

### KMAfinder thresholds
ecoli_kma_threshold = {
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

cdiff_kma_threshold = {
    "tcdA": [90,90,10],
    "tcdB": [90,90,10],
    "tcdC": [90,90,10],
    "cdtAB": [90,90,10],
    "other": [98, 98,10]
}

### Cdiff deletion thresholds
deletion_gt_thresholds = {
    "del117_1": [0.85,1,1], # [IMF, IDV, DP]
    "del330_347_18": [0.01,1,1],
    "del330_347_36": [0.01,1,1],
    "del341_379_39": [0.01,1,1],
    "del313_366_54": [0.01,1,1]
}

### Cdiff deletion thresholds
# [percent of 'N' in region for ambiguous deletions]
deletion_consensus_thresholds = { 
    "del117_1": [0.90], 
    "del330_347_18": [0.90],
    "del330_347_36": [0.80],
    "del341_379_39": [0.80],
    "del313_366_54": [0.80]
}

### AMRfinderplus thresholds

ecoli_amr_threshold = {
    "tet": [90, 95],
    "bla": [90, 95],
    "aac": [90, 95],
    "sul": [90, 95],
    "erm": [90, 95],
    "cat": [90, 95],
    "qnr": [90, 95],
    "mph": [90, 95],
    "van": [90, 95],
    "aad": [90, 95],
    "aph": [90, 95],
    "other": [90, 95]  # default fallback
}

cdiff_amr_threshold = {
    "tet": [90, 95],
    "other": [90, 95]  # default fallback
}


def get_threshold(template_name: str, thresholds: Dict[str, List[int]]) -> List[int]:
    """
    Returns the coverage and identity threshold for a given gene.

    Args:
        template_name (str): Name of the template (gene) from the .res file.
        thresholds (Dict[str, List[int]]): Dictionary of gene thresholds.

    Returns:
        List[int]: A list of two integers: [coverage_threshold, identity_threshold].
    """
    for key in thresholds:
        if key in template_name:
            return thresholds[key]
    return thresholds["other"]

def get_kma_thresholds_for_species(organism_name: str) -> Dict[str, List[int]]:
    if organism_name.strip().lower() in ["e.coli", "e coli", "escherichia coli"]:
        return ecoli_kma_threshold
    elif organism_name.strip().lower() in ["c.diff", "c diff", "clostridium difficile", "clostridioides difficile"]:
        return cdiff_kma_threshold
    raise ValueError(f"No KMA thresholds for: {organism_name}")

def get_deletion_threshold(deletion_key: str, thresholds: Dict[str, List[float]]) -> List[float]:
    """
    Returns the [IMF, IDV, DP] threshold list for a given deletion ID.
    
    Args:
        deletion_key (str): Deletion ID, e.g. 'del330_347_18'

    Returns:
        List[float]: List containing [IMF, IDV, DP] thresholds

    Raises:
        ValueError if the key is not found
    """
    for key in thresholds:
        if key in deletion_key:
            return thresholds[key]
    raise ValueError(f"No deletion thresholds found for: {deletion_key}")


def get_amr_thresholds_for_species(organism_name: str) -> Dict[str, List[int]]:
    if organism_name.strip().lower() in ["e.coli", "e coli", "escherichia coli"]:
        return ecoli_amr_threshold
    elif organism_name.strip().lower() in ["c.diff", "c diff", "clostridium difficile", "clostridioides difficile"]:
        return cdiff_amr_threshold
    raise ValueError(f"No AMRFinder thresholds for: {organism_name}")