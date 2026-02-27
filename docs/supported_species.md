# Pipeline Supported Species

<details>
<summary><strong>Click to expand: Analysis supported </strong></summary>

```
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
    assemblers: shovill

kmeraligner:
    database: elmDB
    options : -ID 80 -1t1 -cge

### Typical of C.Difficile
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

# Typical of E.Faecalis and E.Faecium
LREfinder:
    database : [elmDB]

# Typical of E.Coli
chtyper:
    database: fumCH_db 

custom_blaster:
    options: "-perc_identity 90.0"
    assemblers: shovill
    database : OXAndm

# Typical of K.Pneumoniae
kleborate:
    options: --preset kpsc
    assemblers: shovill  

#  Typical of N.Meningitis
meningotype:
    assemblers: shovill

# Typical of S.Aureus
spatyper:
  assemblers: [spades]    

# S.Enterica
seqsero2: true

sistr:
  assemblers: shovill  

``` 

</details>

## E.Coli
Bla bla
### Genes A & B

## C_difficile  

### Complex workflow etc..

## E_Faecalis & E_Faecium  

## K_pneumoniae  

## N_meningitidis  

## S_aureus  

## S_enterica


