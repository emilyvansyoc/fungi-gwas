## Mendelian randomization analysis using MRSample
# EVS 5/7/2024, updated 3/2025 for updated FAVs

#library(remotes)
library(tidyverse)
#install_github("MRCIEU/TwoSampleMR")
#remotes::install_github('MRCIEU/ieugwasr')
library("TwoSampleMR")
library("ieugwasr")

## the built-in server access is broken; this must be fed manually for each command
# details: https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication
mytoken <- "mytoken"

#### ----- get exposures ----

# need exposures in a specific format
#load("R/gwas-output/sigs_SNPandStructural.RData")
# load 2024 FAVs
load("data/Nov2024_sigs_SNPandStructural.RData")
newsigs <- newsigs %>% mutate(idtaxa = paste(ID, taxa, sep = ":"))
favs <- newsigs %>% distinct(idtaxa)
# need standard errors
sigtaxa <- unique(newsigs$taxa)
# not shared on Github (huge)
load("restricted/effectsizes_allSNP_allvariants_alltaxa.RData")

### format
sigs <- alleffects %>% 
  mutate(taxa = str_remove_all(taxa, "allsnp_allvariants\\."),
         taxa = str_remove_all(taxa, "\\.glm\\.linear")) %>% 
  mutate(idtaxa = paste(ID, taxa, sep = ":")) %>% 
  filter(idtaxa %in% favs$idtaxa)

# remove giant object
rm(alleffects)

# for now just get the heart disease sigs
#hd <- sigs %>% filter(ID %in% c("rs12149890", "rs12929586"))

# format for package
instr <- TwoSampleMR::format_data(sigs, 
                                  snp_col = "ID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "A1_FREQ",
                                  effect_allele_col = "ALT",
                                  other_allele_col = "REF",
                                  pval_col = "P",
                                  phenotype = "taxa"
) # duplicated SNP from Candida and Saccharomyces

# check to make sure all instruemtns are formatted correctly
unique(instr$mr_keep.exposure)

# check for our FAVs in the LD lookup panel
ldref <- ld_reflookup(rsid = instr$SNP, opengwas_jwt = mytoken)
nomatch <- setdiff(unique(sigs$ID), ldref) # 15 FAVs are not in the LD clumping panel
sigs %>% filter(ID %in% nomatch) # some of these are Kaz FAVs -> but are they in LD with other Kaz FAVs?
# all Kaz nomatch are on chr 16 and are structural
# all Capnodiales nomatch are on chr 5 - check for other FAVs on chr 5 around 152600000
# one Capnodiales structural on chrom 10
# two Pleosporales structural on Chr 1

# check LD 
sigs %>% group_by(taxa) %>% count(X.CHROM) # 3 variants on chromosome 10 and 15 on chromosome 1 - OK to lose one
sigs %>% filter(X.CHROM == 5) # tons
sigs %>% filter(X.CHROM == 16) # tons



## get independent variants using their clumping algorithm
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
ieugwasr::api_status()
instr.cl <- clump_data(instr)

# try using IEUGWAS clumping instead
instr1 <- instr %>% rename(rsid = SNP, pval = pval.exposure, id = exposure)
# run this:
options('ieugwas.api' = "https://api.opengwas.io/api/", 'ieugwasr_api' = "https://api.opengwas.io/api/")
instr.cl <- ld_clump(instr1, opengwas_jwt = mytoken) # runs if options are set properly
# run with looser thresholds
instr.cl.loose <- ld_clump(instr1, clump_r2 = 0.00001, clump_kb = 1000, opengwas_jwt = mytoken) # still only 11

## ---- outcomes ----

# get outcome data - runs for a while
ao <- TwoSampleMR::available_outcomes(opengwas_jwt = mytoken)

## get outcomes related to ischemia and coronary heart disease
cvd <- ao %>% 
  filter(str_detect(trait, "schemi"))

myoutcomes <- c("Major coronary heart disease event", "Emergency coronary revascularization (for ACS) (no controls excluded)", "Coronary artery disease (SPA correction)", "Coronary artery disease (Firth correction)", "Major coronary heart disease event excluding revascularizations", "Coronary heart disease", "Coronary atherosclerosis (no controls excluded)", "Coronary artery disease", "Ischemic Heart Disease", "Ischemic heart diseases")

# subset
cvd <- ao %>% 
  filter(trait %in% myoutcomes)

# get outcome data for the SNPs and phenotypes of interest
outc <- extract_outcome_data(
  #snps = instr.clump$SNP,
  snps = instr.cl$rsid,
  outcomes = unique(cvd$id),
  opengwas_jwt = mytoken
)
## some of these have very few SNPs associated with them; remove
outc <- outc %>% 
  filter(!outcome %in% c("Coronary heart disease || id:ieu-a-9"))

# harmonize data
hdat <- harmonise_data(exposure_dat = instr.cl %>% rename(SNP = rsid, exposure = id),
                       outcome_dat = outc, action = 1)

# sanity check
hdat %>% filter(mr_keep == FALSE)

# prune multiple outcomes by choosing the exposure-outcome combination with the highest expected power
# (largest sample size for the outcome)
#pdat <- power_prune(dat = wncase, method = 1, dist.outcome = "binary")
# I don't think this is applicable to us because none of the studies have fungi as an exposure 
# and the code isn't working

# need to prune by hand according to github issue: https://github.com/MRCIEU/TwoSampleMR/issues/80
# this has been a bug since 2019....
# because they don't include `ncase` in extract_outcome_data
wncase <- hdat %>% 
  left_join(ao %>% dplyr::select(id, trait, ncase, ncontrol), by = c("id.outcome" = "id", "originalname.outcome" = "trait")) 

wncase %>% group_by(outcome, ncase) %>% count() %>% arrange(desc(ncase))
## want to select: 
tokeep <- c("Coronary artery disease || id:ebi-a-GCST005195",
            "Ischemic heart diseases || id:finn-b-I9_ISCHHEART",
            "Coronary atherosclerosis (no controls excluded) || id:finn-b-I9_CORATHER_EXNONE",
            "Coronary heart disease || id:ebi-a-GCST000998")

# prune by hand
#mydat <- outc %>% filter(outcome %in% tokeep)
# retain all >20,000 and all 11 instruemtns
tokeep <- wncase %>% group_by(outcome, ncase) %>% count() %>% arrange(desc(ncase)) %>% filter(ncase > 20000 & n == 11) %>% ungroup %>% distinct(outcome) # 8
# filter
hdat <- harmonise_data(exposure_dat = instr.cl %>% rename(SNP = rsid, exposure = id),
                       outcome_dat = outc %>% filter(outcome %in% tokeep$outcome))

# select by hand to remove redudances
#mydat <- hdat %>% 
# filter(!outcome %in% c(""))

# run MR; unclear why one exposure is left out when it is run all in one command; run each exposure separately 
exps <- unique(hdat$exposure)
outres <- data.frame()
for(i in 1:length(exps)) {
  sub <- hdat %>% filter(exposure == exps[i])
  ressub <- mr(sub)
  outres <- rbind(outres, ressub)
}


#res <- mr(hdat, method_list = "mr_wald_ratio")

# get sigs
sigp <- outres %>% filter(pval < 0.05) # 12; all with Wald ratio; 5 taxa and 6 outcome GWAS


### multiple comparisons; there are 8 outcomes
#bonf.10 = 0.05/8
#sigpadj <- res %>% mutate(padjB = if_else(pval < bonf.10, "sig", "NS"))
sigpadj <- outres %>% 
  group_by(exposure) %>% 
  mutate(padj = p.adjust(pval, method = "bonferroni")) %>% 
  ungroup()
# do any meet that threshold?
#sigpadj %>% filter(padjB == "sig") # two with Kazachstania and one Saccharomycetaceae family
sigpadj %>% filter(padj < 0.05) # two with Kaz and one with Saccharomcyetaeae family
allsigs <- sigpadj %>% filter(padj < 0.05)
# what were the outcome GWAS?
cvd %>% filter(id %in% allsigs$id.outcome)

### adjust for ALL comparisons
adj1 <- outres %>% 
  mutate(padjall = p.adjust(pval, method = "bonferroni")) %>% 
  filter(padjall < 0.05)
# what is the outcome GWAS?

# sensitivity: test heterogeneity
het <- mr_heterogeneity(hdat) # only 3 have enough SNPs available (Capnodiales) - none are significant
het %>% filter(Q_pval < 0.05) 

# horizontal pleiotropy (the SNP influences the outcome through some other pathway than the exposure) 
ple <- mr_pleiotropy_test(hdat) # the 3 with enough SNps (Capnodiales) have NA values
ple %>% filter(pval < 0.05) # not enough SNPs available

# do a test on each SNP
sing <- mr_singlesnp(hdat)
singlesig <- sing %>% filter(p < 0.05) %>% group_by(SNP) %>% count() # four unique SNPs
singlesig <- sing %>% drop_na(p) %>% mutate(padj = p.adjust(p, method = "bonferroni")) %>% filter(padj < 0.05) # still the same SNP associated with Kaz and coronary artery disease 

# leave one out to determine if one SNP is driving the association
loo <- mr_leaveoneout(hdat) # same issue with not enough SNPs except for Capnodiales 
loo %>% filter(p < 0.05) %>% group_by(SNP) %>% count() # none

## make a scatter plot
#mr_scatter_plot(res %>% filter(outcome == "Coronary atherosclerosis (no controls excluded) || id:finn-b-I9_CORATHER_EXNONE"), hdat %>% filter(outcome == "Coronary atherosclerosis (no controls excluded) || id:finn-b-I9_CORATHER_EXNONE"))
mr_scatter_plot(res, hdat)
mr_scatter_plot(res %>% filter(exposure == "Genus_Kazachstania"), hdat %>% filter(exposure == "Genus_Kazachstania"))
# makes a scatter plot for each outcome


## make a forest plot
# the black point is the log odds ratio for outcome per standard deviation increase in exposure
mr_forest_plot(mr_singlesnp(hdat %>% filter(outcome == "Coronary atherosclerosis (no controls excluded) || id:finn-b-I9_CORATHER_EXNONE"))) # not all SNPs were found in this outcome 
mr_forest_plot(sing %>% filter(exposure == "Genus_Kazachstania")) # not enough SNPs

# make a funnel plot to assess reliability of the MR analysis
mr_funnel_plot(singlesig) # not enough SNPs

# (don't run): run all analyses, sensitivity analyses, and plots and makes an html report in the working directory
# mr_report(hdat)

# test the directionality of the effect
lrat <- generate_odds_ratios(res %>% filter(id.outcome == "ebi-a-GCST005195" & exposure == "Genus_Kazachstania"))
# weighted median, single mode, and weighted mode are all outside of zero
# effect is negative (CI's are < 0, OR are < 1, 'b' statistic is negative)

# get odds ratios

#get_r_from_lor(lor = lrat$or, af = hdat$eaf.outcome, ncase = hdat$samplesize.outcome, prevalence = )
# try another way to get r2
hkaz <- hdat %>% filter(exposure == "Genus_Kazachstania") %>% filter(mr_keep == TRUE)
r2 <- get_r_from_bsen(b = hkaz$beta.outcome, se = hkaz$se.outcome, n = hkaz$samplesize.outcome)
#hdatr <- hdat %>% cbind(r2) %>% dplyr::rename(r.outcome = r2)%>% filter(outcome == "Coronary atherosclerosis (no controls excluded) || id:finn-b-I9_CORATHER_EXNONE")  %>% mutate(samplesize.exposure = 125)

#directionality_test(hkaz)


# get sigs
sig1 <- res1 %>% filter(pval < 0.05) # 6 are significant; check similarity of outcomes
finn.r2 <- get_r_from_lor(lor = 0.97, af = 0.276, ncase = 23363, ncontrol = 195429, prevalence = 23363/(23363+195429))
finn.r2e <- get_r_from_bsen(b=1.115, se = 0.21, n=125)
hdat1$r.outcome = finn.r2
hdat1$r.exposure = finn.r2e
directionality_test(hdat1 %>% filter(id.outcome == "finn-b-I9_CORATHER_EXNONE"))
mr_steiger(p_exp = 5.21e-7,
           p_out = 6.0e-6,
           n_exp = 125,
           n_out = 122733+424528,
           r_exp = finn.r2e,
           r_out = finn.r2)
