from importlib import resources
from pathlib import Path

# Determining system paths
PKG_DIR = resources.files("mmaseq")
DATA_DIR = PKG_DIR / "data"
READ_DIR = DATA_DIR / "reads"
URL_FILE = READ_DIR / "reads.urls"
SHEET_FULL = DATA_DIR / "samplesheet.tsv"
SHEET_SMLL = DATA_DIR / "samplesheet_small.tsv"
PKG_CONFIGS = PKG_DIR / "config"
SPE_CONFIGS = PKG_CONFIGS / "species_configs"
CATALOGUE_PATH = PKG_CONFIGS / "results_catalogue.yaml"
SCREENING_DIR = PKG_CONFIGS / "target_screening"
WORKFLOW_DIR = PKG_DIR / "workflow"
SNAKEFILE =  WORKFLOW_DIR / "Snakefile"
SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
ENVS_DIR = WORKFLOW_DIR / "envs"
CWD = Path.cwd()
