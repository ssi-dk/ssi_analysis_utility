# Mixed Microbial Analysis on sequencing data - MMAseq
[![Snakemake](https://img.shields.io/badge/snakemake-v7.32-blue?logo=snakemake)](https://snakemake.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

MMAseq is a modular Snakemake-driven workflow that utilizes public tools and custom scripts to coordinate and execute analysis if whole genome microbial sequencing data. It is designed with a configuration based architecture in mind, to facilitate control over the flow of analysis. The modular structure and configuration based architecture enables the user to execute species-specific analysis of a wealth of different species from a single execution.

The pipeline utilizes raw reads and assemblies. In case assemblies are missing, the pipeline offers a few options for de-novo assemblies.

## Tools
MMAseq currently includes 30+ different bioinformatic tools, which will expand in the future.
The included tools enables the pipeline to support ***Assembly***, ***characterization through CGE finders***, ***characterization through NCBI tools***, ***Configurable custom SNP and deletion identification***

## Documentation

A complete documentation is available [MMAseq Documentation](https://ssi-dk.github.io/MMASeq/). The documentation includes a full guide to the workflow from installation, tutorials, bacterial species supported etc...

## Contributing
Please to report bugs or enhancement use the Issue Tracker. 

# Citation
When using SSI analysis in published work, please cite the following...


For citations of included algorithms and sub-modules please see the (placeholerreferences).
