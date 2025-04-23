# fungi-gwas
Code for analyses and figures for Van Syoc et al, under review.  

NOTE that most data cannot be shared due to genomics data sharing restrictions. Throughout each analysis script, restricted data is marked with the dummy pathfile "/restricted/". Otherwise, most scripts are fully reproducible and executed by cloning the repository to a local machine and running the scripts. (With the caveat that conda environment YAML files are not currently shared; a conda environment can be created with the softwares installed with default mechanisms. No special software versions were used in these analyses)

## Analyses 

### Amplicon processing 

A series of scripts pulls amplicon data from Sequence Read Archive and processes it to generate an OTU table and taxonomy. First, `get-runinfo.sh` and `fastq-dump.sh` pull data from SRA. Then, `ITS_qualitycontrol.R` runs basic quality control, then `vsearch_ITS.sh` runs OTU generation.

### GWAS (PLINK)

First, fungal data is parsed to prep for GWAS. This is in `ITS-vsearch-to-phyloseq_prepGWAS.R`.

The script `plink-code.sh` runs the basic steps of GWAS using Human Microbiome Data accessed from dbGaP. After preprocessing, the GLM model is run with a pheno file that has one column for the relative abundance of each fungal taxa, and a covariate file with one column for each covariate. PLINK iterates through each column of the pheno file and writes a separate text file for the results of each fungal taxa. 

### post-GWAS 

1. `get_FAVs.R` loads an RData object of all combined GWAS results (giant!) and gets FAVs at various significance levels  
2. `SNPNexus.R` parses the results from SNP annotation done through SNPNexus web browser  
3. `GTEx.R` fetches FAV-eQTLs using the lookup table from GTEx v10  
4. `phewas.R` parses the PheWEB results from a text file (FAVs were searched in the PheWEB database and top hits pulled by hand - strongly do not recommend but no code-able options) and performs a validation analysis using the ieugwasr R package  