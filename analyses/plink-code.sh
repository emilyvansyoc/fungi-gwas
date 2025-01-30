### PLINK code for running mycobiome GWAS 
# EVS 1/2024

# assuming that DIR=/path/to/files and PLINK and PLINK2 are /paths/to/plink/executables
# this reads in the VCF directly from dbGaP files -> no additional preprocessing steps

#### ---- QUALITY CONTROL ----

# read in data from HMP VCF
$PLINK --vcf $VCF --vcf-filter --double-id --make-bed

# make a missing report (individual and SNP-wise)
$PLINK --bfile plink --missing

# remove SNPs with missing rate > 5%
$PLINK --bfile plink --geno 0.05 --make-bed --out snpmiss05

# remove individuals with missing > 10%
$PLINK --bfile snpmiss05 --mind 0.1 --make-bed --out indmiss10

# check sex; use annotated FAM file with sex from Kolde 2018
$PLINK --bfile indmiss10 --check-sex --fam indmiss.fam

# get only autosomes
$PLINK --bfile indmiss10 --chr 1-22 --make-bed --out autosomes05

# remove SNPs with minor allele frequency < 10%
$PLINK --bfile autosomes05 --maf 0.1 --make-bed --out maf10

# remove SNPs outside of hardy-weinburg equilibrium (using default recommended P value)
$PLINK --bfile maf10 --hwe 1e-10 --make-bed --out hwe_maf10

# remove individuals with heterozygosity rates > 3*standard deviation (suggests inbreeding or sample contamination)
## to do this, first rename SNPs with missing IDs (very few but throws off the way PLINK indexes variants to keep/remove)
## then get independent SNPs
$PLINK --bfile hwe_maf10 --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out renamed_hwemaf10
$PLINK --bfile renamed_hwemaf10 --indep-pairwise 50 5 0.2 --out indepSNP
$PLINK --bfile renamed_hwemaf10 --extract indepSNP.prune.in --make-bed --out indepSNP
$PLINK --bfile indepSNP --het --out R_check
# this next step uses an R script from Marees 2018 tutorial to flag individuals with high heterozygosity
Rscript check_heterozygosity_rate.R
Rscript heterozygosity_outliers_list.R
# then use the list of individuals to remove them
cat fail-het-qc.txt | cut -f 1,2 > rem_fail-het-qc.txt

### from this point forward, have TWO working PLINK files: independent SNPs, and all SNPs
# this also adds back in sexes from a fam file
$PLINK --bfile hwe_maf10 --remove rem_fail-het-qc.txt --fam sexes.fam --make-bed --out rem_het3x
$PLINK --bfile indepSNP --remove rem_fail-het-qc.txt --fam sexes.fam --make-bed --out indepSNP_remhet

# finally, remove the hapmap controls
$PLINK --bfile rem_het3x --remove controlsamp.txt --make-bed --out final_remhet3x
$PLINK --bfile indepSNP_remhet --remove controlsamp.txt --make-bed --out final_indepSNP

## ---- POPULATION STRATIFICATION ----

# do PCA in PLINK2
# this is done on only independent SNPs
$PLINK2 --bfile final_indepSNP --pca --out pca_indepSNP_rem-het
# the output file ".eigenvec" is used as covariates

## ---- RUN GLMs ----

### run on all SNPs
# first subset to just the data that we have fungi for (determined in R)
$PLINK --bfile final_remhet3x --extract ids-keep.txt --make-bed --out final_remhet_sub

# this is run in a bash file; basic structure of the GLM for both comparisons
$PLINK2 --bfile $BFILE --glm sex hide-covar --covar $COVAR --pheno $PHENO --num_threads 10 --out $OUT --adjust
# where PHENO is phentest_CLR_rarefied_cortax.txt and COVAR is covars_remhat_cat_nosex_noseqd.txt

