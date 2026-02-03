# Getting Started

## Installation Guide

The latest version of **SSI Analysis Utility** can be installed using [Conda](#conda-installation), [Pip](#pip), or by installing it [from source](#install-from-source).

!!! info

    SSI Analysis Utility currently supports **Linux** and **macOS** only.

---

## Conda Installation

Because the pipeline is developed in **Snakemake**, and Snakemake does **not** fully support micromamba, installation via **Conda** is recommended.  
To ensure a clean and isolated setup, create a fresh conda environment:

```bash
conda create -n ssi_analysis_utility -c conda-forge -c bioconda ssi_analysis_utility
```

Activate the enviroment

```bash
conda activate ssi_analysis_utility
```

### <a name="pip"></a> pip

You can install the latest stable release of SSI Analysis Utility via pip:

```bash
pip install ssi_analysis_utility 
```

### <a name="Install from source"></a> Install from source  

This section is for user who wants to install **SSI Analysis utility** directly from source (e.g., for reproducibility or to debug a release).
However, the user must ensure to have the **Snakemake >=9.1.1** dependency is installed in the environment. 

```bash
conda create -n ssi_analysis_utility -c conda-forge -c bioconda "snakemake>=9.12.0"
git clone https://github.com/your-org/ssi_analysis_utility.git
```

If the user wants to work on a development branch:

```bash
conda create -n ssi_analysis_utility -c conda-forge -c bioconda "snakemake>=9.12.0"
git clone https://github.com/your-org/ssi_analysis_utility.git
cd ssi_analysis_utility
git switch dev    # work on the develop branch
git checkout -b new-feature-branch  # Create a new branch for new changes
``` 