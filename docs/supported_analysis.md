# Supported Analysis

The pipeline is designed to be flexible and modular, allowing users to customize the analysis based on their specific research questions and the characteristics of their data. The workflow includes a range of analysis options that can be applied to both raw reads and assembled genomes, depending on the user's preferences and the requirements of their study. The specific analyses that can be performed are determined by the configuration settings provided in the samplesheet and the species-specific configuration files. These settings allow users to select which tools and databases to use for each analysis, as well as to specify any additional parameters or options that may be relevant for their particular dataset or research question.

## Generic Analysis
The pipeline supports a variety of generic analyses that can be applied to a wide range of bacterial species. Following are listed the analysis options that can be applied to all species, depending on the configuration settings:


```yaml
mlst:
    assemblers: shovill  
    
resfinder:
    options: --species 'Other'

plasmidfinder:
    options: ""

virulencefinder:
    options: ""

serotypefinder:
    options: ""
          
amrfinder:
    assemblers: shovill
    options: ""
    
meningotype:
    assemblers: [shovill, spades, skesa]
```

The options field can be used to specify any additional parameters or options that may be relevant for the analysis, such as the species being analyzed or specific thresholds for gene detection and can be specified following the rules of the respective tool. The assemblers field allows users to specify which assembly tools to use for the analysis, which can be particularly important for certain analyses that may require specific assembly methods to achieve optimal results (in the current states only shovill, spades and skesa are supported). By providing these configuration options, the pipeline allows users to tailor their analysis to their specific needs and research questions, while still maintaining a consistent and standardized workflow across different samples and species.


## Species Specific analysis

In addition to the generic analyses that can be applied to all species, the pipeline also supports species-specific analyses that are tailored to the unique characteristics and research questions associated with specific bacterial species. These analyses are defined in the species-specific configuration files and can include additional tools, databases, and parameters that are relevant for the particular species being studied. By providing support for both generic and species-specific analyses, the pipeline allows users to conduct comprehensive and customized analyses that are optimized for their specific research needs and the characteristics of their data. The following sections provide examples of species-specific analyses that can be performed for different bacterial species, demonstrating the flexibility and versatility of the pipeline in accommodating a wide range of research questions and genomic contexts.

## E.Coli

``` yaml

kma_filter:
    options: --organism 'Escherichia coli'
    database: ecoligenes

chtyper:
    database: fumCH_db 

custom_blaster:
    options: "-perc_identity 90.0"
    assemblers: shovill
    database : OXAndm
```


## C.Difficile  

```yaml
kma_filter:
    options: --organism 'Clostridioides difficile'
    database: Cdiff_Toxins

snp_identifier:
    options: --organism 'Clostridioides difficile' --snp_region_buffer 5
    database: Cdiff_Toxins

deletion_identifier:
    options: --organism 'Clostridioides difficile' --deletion_region_buffer 5 --overlap_fraction 0.5
    database: Cdiff_Toxins
    assemblers: [spades]        # <--- NEW: matches the assembler wildcard

assembly_minimap2:
    options: -ax asm5 --cs
    database: Cdiff_Toxins
    assemblers: [spades]        # used to configure which assembler(s) to align

cdiff_repeat_identifier:
    repeats: TR6 TR10
    combos: TRST
    assemblers: [spades, skesa]

samtools:
    view_options: -q 20
    sort_options: -n

bcftools:
    view_options: -v indels -i 'INFO/DP>10'
```


## E.Faecalis - E.Faecium  

```yaml
LREfinder:
    database : [elmDB]
```

## K.Pneumoniae  
```yaml
kleborate:
    options: --preset kpsc
    assemblers: shovill  
```

## N_meningitidis  

```yaml
meningotype:
    assemblers: shovill
```

## S.Aureus  

# Typical of S.Aureus

```yaml
spatyper:
  assemblers: [spades]   
```
## S.Enterica

```yaml
seqsero2: true

sistr:
  assemblers: shovill  
```


