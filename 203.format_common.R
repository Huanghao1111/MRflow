#1.read data
#2.format data
#3.screen pval
#4.clump data
#5.Save Data at Each Step

dat <- data.table::fread("data/01origin_data/data.gz")
dat <- format_data(dat, id, SNP_col = "rsid",
                   chr_col = "chromosome",pos_col = "base_pair_location",
                   ea_col = "effect_allele", nea_col = "other_allele",
                   eaf_col = "effect_allele_frequency",
                   beta_col = "beta", se_col = "standard_error",
                   pval_col = "p_value", n_col="samplesize")
dat$Phenotype <- "Phenotype"
data.table::fwrite(dat, "data/02format_data/data.gz")
#
dat <- subset(dat, pval<5e-8)
data.table::fwrite(dat, "data/03filter_data/data.gz")
dat$pval.exposure <- dat$pval
#clump data
dat <- TwoSampleMR::clump_data(dat, clump_p1 = 5e-8, bfile = "LD/bfile/EUR", plink_bin="LD/plink_win64/plink.exe")
dat$id.exposure = NULL
dat$pval.exposure = NULL
data.table::fwrite(dat, "data/04clump_data/data.gz")
