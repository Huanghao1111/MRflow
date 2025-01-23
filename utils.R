library(dplyr)
#### format data ####
format_data <- function(dat, id, SNP_col="rsid", chr_col="chr", pos_col="position",
                        ea_col="ea", nea_col="nea", eaf_col="eaf",
                        beta_col="beta", se_col="se",
                        pval_col="p",
                        n_col="n"
)
{
  columns_data <- dplyr::select(dat, SNP=SNP_col, chr=chr_col, pos=pos_col,
                                effect_allele=ea_col, other_allele=nea_col, eaf=eaf_col,
                                beta=beta_col,se=se_col,
                                pval=pval_col, samplesize=n_col)

  columns_data$id <- id
  #remove na data
  columns_data <- subset(columns_data, SNP != "")
  columns_data <- columns_data[!grepl(":", columns_data$SNP), ]
  columns_data <- columns_data[grepl("rs", columns_data$SNP), ]
  if(!is.numeric(columns_data$p)){
    columns_data$pval <- as.numeric(columns_data$pval)
  }
  #change X/Y to 22/23
  if(!is.numeric(columns_data$chr)){
    columns_data$chr <- as.numeric( ifelse(columns_data$chr == "X" | columns_data$chr=="Y",
                                           ifelse(columns_data$chr == "X", "22", "23"),
                                           columns_data$chr
    ))
  }
  return(columns_data)
}


#### get_gcst_msg ####
get_gcst_msg<- function(file_path, gcst_id){
  if(grepl(".xlsx", file_path)){
    a=openxlsx::read.xlsx(file_path)
  }else{
    a = data.table::fread(file_path)
  }
  return(subset(a, accessionId==gcst_id))
}
#### F ####
cal_F <- function(dat){
  dat$R2 <- 2 * (1-dat$eaf) * dat$eaf * dat$beta^2
  dat$F <- (dat$samplesize-2)*(dat$R2 / (1-dat$R2))
  dat
}

#### fdr  ####
cal_fdr <- function(mr_results){

  mr_results <- mr_results %>%
    mutate(original_order = row_number())


  mr_results <- mr_results %>%
    group_by(method) %>%
    mutate(p_fdr = p.adjust(pval, method = "fdr")) %>%
    ungroup()


  mr_results <- mr_results %>%
    arrange(original_order) %>%
    select(-original_order) # 移除临时列

  return(mr_results)
}
