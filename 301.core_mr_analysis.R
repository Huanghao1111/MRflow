source("utils.R")
#' @title Two-Sample Mendelian Randomization Analysis Pipeline
#' @description This function performs a two-sample Mendelian Randomization (MR) analysis using specified exposure and outcome datasets. It is designed to streamline the MR analysis process, facilitating the investigation of causal relationships between genetic variants and traits.
#' @param exposure_path A string specifying the file path to the exposure data. This file should contain relevant data for the exposure variable(s) used in the analysis.
#' @param outcome_path A string specifying the file path to the outcome data. This file should contain relevant data for the outcome variable(s) used in the analysis.
#' @param output_path A string specifying the directory path where the analysis results will be saved. This should be a writable location.
tsmr_analysis <- function(exposure_path, outcome_path, output_path){
  #exposure data
  exposure_data <- data.table::fread(exposure_path)
  exposure_data <- TwoSampleMR::format_data(data.frame(exposure_dat), type="exposure")
  #outcome data
  outcome_data <- data.table::fread(outcome_path)
  outcome_data <- TwoSampleMR::format_data(data.frame(outcome_data),
                                           snps=exposure_data$SNP, type="outcome")
  #harmonise data
  harmonise_data <- TwoSampleMR::harmonise_data(exposure_data, outcome_data, action = 2)
  #mr analysis
  mr_results <- TwoSampleMR::mr(harmonise_data, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  #or
  mr_results <- TwoSampleMR::generate_odds_ratios(mr_results)
  #fdr
  mr_results <- cal_fdr(mr_results)
  #Heterogeneity Test
  heterogeneity_test = TwoSampleMR::mr_heterogeneity(harmonise_data)
  #Pleiotropy test
  pleiotropy_test = TwoSampleMR::mr_pleiotropy_test(harmonise_data)
  #save results
  if(!dir.exists(output_path)){
    dir.create(output_path)
  }
  data.table::fwrite(mr_results, file.path(output_path, "mr_results.csv"))
  data.table::fwrite(heterogeneity_test, file.path(output_path, "heterogeneity_test.csv"))
  data.table::fwrite(pleiotropy_test, file.path(output_path, "pleiotropy_test.csv"))
}
#Forward
tsmr_analysis("data/04clump_data/1400met_all.gz", "data/02format_data/ptsd/PTSD.gz", "results/1400_ptsd/")
tsmr_analysis("data/04clump_data/338met_all.gz", "data/02format_data/ptsd/PTSD.gz", "results/338_ptsd/")
tsmr_analysis("data/04clump_data/233met_all.gz", "data/02format_data/ptsd/PTSD.gz", "results/233_ptsd/")

tsmr_analysis("data/04clump_data/1400met_all.gz", "data/02format_data/ptsd/finngen_R11_F5_PTSD.gz", "results/1400_ptsd_finn/")
tsmr_analysis("data/04clump_data/338met_all.gz", "data/02format_data/ptsd/finngen_R11_F5_PTSD.gz", "results/1400_ptsd_finn/")
tsmr_analysis("data/04clump_data/233met_all.gz", "data/02format_data/ptsd/finngen_R11_F5_PTSD.gz", "results/1400_ptsd_finn/")

#reverse
tsmr_analysis("data/04format_data/PTSD_clump.gz", "data/02format_data/233met/GCST90302016.gz", "results/ptsd_GCST90302016/")
tsmr_analysis("data/04format_data/PTSD_clump.gz", "data/02format_data/233met/GCST90302030.gz", "results/ptsd_GCST90302030/")
tsmr_analysis("data/04format_data/PTSD_clump.gz", "data/02format_data/233met/GCST90302097.gz", "results/ptsd_GCST90302097/")
tsmr_analysis("data/04format_data/PTSD_clump.gz", "data/02format_data/233met/GCST90302150.gz", "results/ptsd_GCST90302150/")
tsmr_analysis("data/04format_data/PTSD_clump.gz", "data/02format_data/233met/GCST90302167.gz", "results/ptsd_GCST90302167/")

tsmr_analysis("data/04format_data/PTSD_clump.gz", "data/02format_data/338met/GCST90026178.gz", "results/ptsd_GCST90026178/")

tsmr_analysis("data/04format_data/PTSD_clump.gz","data/02format_data/1400met/GCST90200371.gz", "results/ptsd_GCST90200371/")
tsmr_analysis("data/04format_data/PTSD_clump.gz","data/02format_data/1400met/GCST90200716.gz", "results/ptsd_GCST90200716/")
tsmr_analysis("data/04format_data/PTSD_clump.gz","data/02format_data/1400met/GCST90200877.gz", "results/ptsd_GCST90200877/")
tsmr_analysis("data/04format_data/PTSD_clump.gz","data/02format_data/1400met/GCST90200559.gz", "results/ptsd_GCST90200559/")
tsmr_analysis("data/04format_data/PTSD_clump.gz","data/02format_data/1400met/GCST90200099.gz", "results/ptsd_GCST90200099/")
tsmr_analysis("data/04format_data/PTSD_clump.gz","data/02format_data/1400met/GCST90200207.gz", "results/ptsd_GCST90200207/")
tsmr_analysis("data/04format_data/PTSD_clump.gz","data/02format_data/1400met/GCST90199905.gz", "results/ptsd_GCST90199905/")
tsmr_analysis("data/04format_data/PTSD_clump.gz","data/02format_data/1400met/GCST90200947.gz", "results/ptsd_GCST90200947/")
tsmr_analysis("data/04format_data/PTSD_clump.gz","data/02format_data/1400met/GCST90199944.gz", "results/ptsd_GCST90199944/")

