source("utils.R")
library(dplyr)
#https://pgc.unc.edu/for-researchers/download-results/
#https://doi.org/10.6084/m9.figshare.26349322
dat <- data.table::fread("data/01origin_data/ptsd/eur_ptsd_pcs_v4_aug3_2021.vcf")
dat$se <- 1 / sqrt((2 * dat$FREQ) * (1 - dat$FREQ) * (dat$NEFF + (dat$Z^2)))
dat$beta <- dat$Z * dat$se
dat <- format_data(dat, "PTSD", SNP_col = "ID",
                   chr_col = "#CHROM",pos_col = "POS",
                   ea_col = "A1", nea_col = "A2",
                   eaf_col = "FREQ",
                   beta_col = "beta", se_col = "se",
                   pval_col = "P", n_col="NEFF")
dat$Phenotype <- "PTSD"
data.table::fwrite(dat, "data/02format_data/ptsd/PTSD.gz")
#
dat <- subset(dat, pval<5e-8)
data.table::fwrite(dat, "data/03filter_data/PTSD.gz")
dat$pval.exposure <- dat$pval
#clump data
dat <- TwoSampleMR::clump_data(dat, clump_p1 = 5e-8, bfile = "LD/bfile/EUR", plink_bin="LD/plink_win64/plink.exe")
dat$id.exposure = NULL
dat$pval.exposure = NULL
data.table::fwrite(dat, "data/04clump_data/PTSD_clump.gz")

#https://www.finngen.fi/en/access_results
#https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv
#https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_F5_PTSD.gz
dat <- data.table::fread("data/01origin_data/ptsd/finngen_R11_F5_PTSD.gz")
dat$samplesize <- 3005+403817
dat <- format_data(dat, "finngen_R11_F5_PTSD", SNP_col = "rsids",
                   chr_col = "#chrom",pos_col = "pos",
                   ea_col = "alt", nea_col = "ref",
                   eaf_col = "af_alt",
                   beta_col = "beta", se_col = "sebeta",
                   pval_col = "pval", n_col="samplesize")
dat$Phenotype <- "Post-traumatic stress disorder"
data.table::fwrite(dat, "data/02format_data/ptsd/finngen_R11_F5_PTSD.gz")

