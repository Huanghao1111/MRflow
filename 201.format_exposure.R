#Read and Format Exposure Data
#This process involves several time-consuming steps:

# 1. Format the Data
#Begin by organizing and structuring the raw data into a consistent format suitable for analysis.

# 2. Screen Data for P-Values
# Filter the data based on p-values to focus on statistically significant results.

# 3. Clump the Data
#Group related data points to reduce redundancy and enhance analysis efficiency.

# 4. Save Data at Each Step
#Ensure that data is saved after each step to prevent data loss and facilitate future reference.


source("utils.R")
library(dplyr)

#### 1400met ####
files <- list.files("data/01origin_data/1400met/", full.names = T)
for(file in files){
  dat <- data.table::fread(file)
  id <- gsub(".h.tsv.gz", "", basename(file))
  msg <- get_gcst_msg("url_export/PMID36635386_studies_export.tsv", id)
  dat$samplesize <- as.numeric(gsub(",", "", gsub(" European ancestry individuals", "", msg$initialSampleDescription)))
  dat <- format_data(dat, id, SNP_col = "rsid",
                    chr_col = "chromosome",pos_col = "base_pair_location",
                    ea_col = "effect_allele", nea_col = "other_allele",
                    eaf_col = "effect_allele_frequency",
                    beta_col = "beta", se_col = "standard_error",
                    pval_col = "p_value", n_col="samplesize")
  dat$Phenotype <- msg$reportedTrait
  data.table::fwrite(dat, file.path("data/02format_data/1400met/", paste0(id, ".gz")))
  #filter data
  dat <- subset(dat, pval<1e-5)
  data.table::fwrite(dat, file.path("data/03filter_data/1400met/", paste0(id, ".gz")))
  dat$pval.exposure <- dat$pval
  #clump data
  dat <- TwoSampleMR::clump_data(dat, clump_p1 = 1e-5, bfile = "LD/bfile/EUR", plink_bin="LD/plink_win64/plink.exe")
  dat$id.exposure = NULL
  dat$pval.exposure = NULL
  data.table::fwrite(dat, file.path("data/04clump_data/1400met/", paste0(id, ".gz")))
}
#merge 1400 clump data for exposure analysis
files = list.files("data/04clump_data/1400met/", full.names = T)
merged_data <- data.frame()
for (file in files) {
  data <- data.table::fread(file)
  merged_data <- rbind(merged_data, data)
}
data.table::fwrite(merged_data, "data/04clump_data/1400met_all.gz")
#### 233met ####
files <- list.files("data/01origin_data/233met/", full.names = T)
for(file in files){
  dat <- data.table::fread(file)
  id <- gsub(".h.tsv.gz", "", basename(file))
  msg <- get_gcst_msg("url_export/PMID38448586_studies_export.tsv", id)
  dat$samplesize <- msg$samplesize
  dat <- format_data(dat, id, SNP_col = "rsid",
                     chr_col = "chromosome",pos_col = "base_pair_location",
                     ea_col = "effect_allele", nea_col = "other_allele",
                     eaf_col = "effect_allele_frequency",
                     beta_col = "beta", se_col = "standard_error",
                     pval_col = "p_value", n_col="samplesize")
  dat$Phenotype <- msg$reportedTrait
  data.table::fwrite(dat, file.path("data/02format_data/233met/", paste0(id, ".gz")))
  #filter data
  dat <- subset(dat, pval<5e-8)
  data.table::fwrite(dat, file.path("data/03filter_data/233met/", paste0(id, ".gz")))
  dat$pval.exposure <- dat$pval
  #clump data
  dat <- TwoSampleMR::clump_data(dat, clump_p1 = 5e-8, bfile = "LD/bfile/EUR", plink_bin="LD/plink_win64/plink.exe")
  dat$id.exposure = NULL
  dat$pval.exposure = NULL
  data.table::fwrite(dat, file.path("data/04clump_data/233met/", paste0(id, ".gz")))
}
#merge 233 clump data for exposure analysis
files = list.files("data/04clump_data/233met/", full.names = T)
merged_data <- data.frame()
for (file in files) {
  data <- data.table::fread(file)
  merged_data <- rbind(merged_data, data)
}
data.table::fwrite(merged_data, "data/04clump_data/233met_all.gz")



#### 338met ####
files <- list.files("data/01origin_data/338met/", full.names = T)
for(file in files){
  dat <- data.table::fread(file)
  id <- unlist(strsplit(gsub(".h.tsv.gz", "", basename(file)), "-"))[2]

  msg <- get_gcst_msg("url_export/PMID33437055_studies_export.tsv", id)
  dat$samplesize <- as.numeric(gsub(",", "", gsub(" European", "", msg$discoverySampleAncestry)))
  dat <- format_data(dat, id, SNP_col = "hm_rsid",
                     chr_col = "hm_chrom",pos_col = "hm_pos",
                     ea_col = "hm_effect_allele", nea_col = "hm_other_allele",
                     eaf_col = "hm_effect_allele_frequency",
                     beta_col = "hm_beta", se_col = "standard_error",
                     pval_col = "p_value", n_col="samplesize")
  dat$Phenotype <- msg$reportedTrait
  data.table::fwrite(dat, file.path("data/02format_data/338met/", paste0(id, ".gz")))
  #filter data
  dat <- subset(dat, pval<5e-5)
  dat$Phenotype <- msg$reportedTrait
  data.table::fwrite(dat, file.path("data/03filter_data/338met/", paste0(id, ".gz")))
  dat$pval.exposure <- dat$pval
  #clump data
  dat <- TwoSampleMR::clump_data(dat, clump_p1 = 5e-5, bfile = "LD/bfile/EUR", plink_bin="LD/plink_win64/plink.exe")
  dat$id.exposure = NULL
  dat$pval.exposure = NULL
  data.table::fwrite(dat, file.path("data/04clump_data/338met/", paste0(id, ".gz")))
}
#merge 338 clump data for exposure analysis
files = list.files("data/04clump_data/338met/", full.names = T)
merged_data <- data.frame()
for (file in files) {
  data <- data.table::fread(file)
  merged_data <- rbind(merged_data, data)
}
data.table::fwrite(merged_data, "data/04clump_data/338met_all.gz")

