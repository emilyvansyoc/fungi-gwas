# fungi-gwas
Code for analyses and figures for Van Syoc et al, under review

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