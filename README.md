# fungi-gwas
Code for analyses and figures for Van Syoc et al, under review.  

[![DOI](https://zenodo.org/badge/923057498.svg)](https://doi.org/10.5281/zenodo.15659049)

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
5. `MR.R` performs Mendelian Randomization using the Two-Sample MR package  
6. `all_Kazachstania.R` takes the phyloseq objects of ITS data from public datasets (see Re-analysis of public ITS data below) and does PERMANOVA on each, then plots (creates Figure S4)  

### Figure generation 
1. `Figure1.R`  
2. `Figure2.R`  
3. `Figure3.R`  (includes minor stats for Kazachstania comparisons) 
4-7. Supplementary Figures: S1-S4

### Re-analysis of public ITS data  

The script `my-its-pipeline.sh` in the folder `ITS_reanalysis` takes a Bioproject accession number and prefix of Biosample runs (i.e., "SRR" for American data, "ERR" for data uploaded to the European Nucleotide Archive that is also hosted on SRA). The pipeline then retrieves the data from SRA, runs primer removal if a fasta file of primer sequences is supplied, and then runs the standard VSEARCH pipeline for quality control, OTU clustering, and taxonomic assignment on a user-supplied database. It finally runs an Rscript to convert VSEARCH output into a phyloseq object and writes an R data file. This relies on a micromamba environment but can be replaced with any conda environment with the necessary softwares, or running without conda assuming that each software is accessible in PATH.  
The other scripts in this folder are called by the main 'my-its-pipeline.sh' script. 
*Note: this is currently set up to use only paired-end data and will fail if single-end data is supplied*
