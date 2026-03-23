# *Clostridioides difficile*

This species-specific analysis pipeline is used to analyze samples belonging to the *Clostridioides difficile* species and determine the presence of predefined canonical sequence features, including repeat-associated signals, SNPs, and canonical deletions such as 18 bp, 39 bp, and 54 bp deletion events in target loci.

The workflow combines alignment-based read evidence, pileup-based indel evidence, consensus-sequence ambiguity, and optional assembly support to classify variants in a biologically interpretable way.

## Requirements
### Database
Several distinct databases are required to determine repeat-associated signals, SNPs, and deletions. 
- Using a customized [python script](https://github.com/ssi-dk/ssi_analysis_utility/blob/dev/workflow/scripts/genbank_fetcher.py), the desired regions of one or multiple loci/locus can be downloaded directly from [genbank](https://www.ncbi.nlm.nih.gov/genbank/) based on their accession numbers. The sequences of interest can be downloaded in several ways, depending on the information provided.
  1. The loci/locus can be downloaded in its entirety by simply providing the accession number(s) within the [metafile](https://github.com/ssi-dk/ssi_analysis_utility/blob/dev/config/Metadata/meta_locus.tsv).
  2. Any coding sequencing within a particular locus can be downloaded based on the locus-specific coordinates, within the [metafile](https://github.com/ssi-dk/ssi_analysis_utility/blob/dev/config/Metadata/meta_cds.tsv)
  3. The loci/locus/genomic region can be downloaded based on the accession number and the chromosomal coordinates within the [metafile](https://github.com/ssi-dk/ssi_analysis_utility/blob/dev/config/Metadata/meta_acc.tsv).
-  The repeat-associated signals are defined within a curated database 
### Bioinformatical tools
- The component is alignment-based using the tool [KMA](https://bitbucket.org/genomicepidemiology/kma).
- The pileup-based evidence is generated from the aligned files ([SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)) and processed using the [htslib](https://www.htslib.org/) library involving a suite of programs, i.e [samtools](https://github.com/samtools/samtools) and [bcftools](https://samtools.github.io/bcftools/)
- The alignment and pileup-based evidence is supported by *de novo* generated assemblies obtained either from [SKESA](https://github.com/ncbi/SKESA), [spades](https://github.com/ablab/spades) or [shovill](https://github.com/tseemann/shovill)

### Snakemake environments
- The alignments creating the Sequence Alignment Map format and the consensus sequence utilize the [KMA](https://github.com/ssi-dk/ssi_analysis_utility/blob/dev/workflow/envs/kmeraligner.yaml)
- The pileup-based functionality for filtering, sorting and indexing of the alignment files, followed by pileup generation, genotype calling and indel filtering are all based on the [htslib](https://github.com/ssi-dk/ssi_analysis_utility/blob/dev/workflow/envs/htslib.yaml) library environment. 
- The *de novo* assembly tools are all based on an environment for [shovill](https://github.com/ssi-dk/ssi_analysis_utility/blob/dev/workflow/envs/shovill.yaml) which includes the SKESA and spades internally.
- The following downstream analysis utilizes customized Python scripts wrangling all of these generated results (described later) and are all based on the same [python environment](https://github.com/ssi-dk/ssi_analysis_utility/blob/dev/workflow/envs/python_functions.yaml)
### Configuration files
To accurately filter the generated alignment files to identify the genes of interest, deletetorioys regions, specific SNP and repeat-associated signals, several configuration files used by the customized Python scripts are required.


## Gene filtering
- The component is alignment-based using the tool [KMA](https://bitbucket.org/genomicepidemiology/kma).




## Tools used for analysis
Each component can be run on each sample individually using one snakemake command, replacing the string passed to the **--config sample_name=" "** with the appropriate dataset name. The provided **component_name=** takes as an argument *<component_name>__<version_number>*. The component name aligns with the GitHub repo name, which is structured like *bifrost_<component_name>* (e.g. *bifrost_sp_ecoli* -> component name *sp_ecoli*), and the version number aligns with the current [GitHub tag](https://github.com/ssi-dk/bifrost_sp_ecoli/tags) / or conda environment [version](https://github.com/ssi-dk/bifrost_sp_ecoli/blob/main/setup.py) (e.g. *v.0.0.2*) defined during the bifrost component setup. 
```bash
insert cmd
```

## Repeat identification

## SNP identification

## Deletion identification

## Results

### KMA
### SPADES
### Repeat Identification
### SNP 
### Spades

One example of the aligned kma results *kma.res*
```bash
#Template       Score   Expected        Template_length Template_Identity       Template_Coverage       Query_Identity  Query_Coverage  Depth   q_value p_value
1__wzx__O157__JH959508     78656            1730            1392          100.00          100.00          100.00          100.00           56.60        73613.56        1.0e-26
2__wzy__O157__JH953200     24223            1532            1185          100.00          100.00          100.00          100.00           20.27        19989.30        1.0e-26
5__fliC__H7__AF228487     359574            1726            1758          100.00          100.00          100.00          100.00          204.77        354427.40       1.0e-26
12__eae__eae-42__AF071034         440253            2546            2805          100.00          100.00          100.00          100.00          157.42        432671.71       1.0e-26
32__ehxA__ehxA-3__AB011549        134647            3574            2997           99.97          100.00           99.97          100.00           45.08        124291.98       1.0e-26
6__stx2__stx2-a__X07865   305905            1280            1241           99.92          100.00           99.92          100.00          246.78        302085.84       1.0e-26
7__stx2__stx2-c__AB015057          70033            1552            1241          100.00          100.00          100.00          100.00           56.33        65510.17        1.0e-26
```
