import json
import os
import subprocess
from pathlib import Path
import shutil 
import pandas as pd
import yaml

conda_bin = shutil.which("conda") 
if not conda_bin:
    raise RuntimeError("conda not found — ensure conda is installed in the current environment")

deploy_dir = snakemake.input.deploy_dir
versions_file = snakemake.output.versions_file

conda_dir = Path(f"{deploy_dir}/conda").absolute()
versions_list = []

for conda_path in conda_dir.iterdir():
    if conda_path.suffix != ".yaml":
        continue
    with open(conda_path) as f:
        try:
            conda_env = yaml.safe_load(f)
        except yaml.YAMLError as e:
            print(e)
            continue

    prefix = conda_path.with_suffix('')
    if not prefix.is_dir():
        continue
    result = subprocess.run(
        [conda_bin, "list", "--prefix", str(prefix), "--json"],
        capture_output=True, text=True, check=True
    )
    for pkg in json.loads(result.stdout):
        if pkg["name"] == conda_env.get("name"):
            versions_list.append({"Tool": pkg["name"], "Version": pkg["version"]})

pd.DataFrame(versions_list).to_csv(versions_file, sep="\t", index=False)
