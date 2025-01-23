#This script can use to build batch url in catalog.

#https://www.ebi.ac.uk/gwas/publications/36635386
export_data <- data.table::fread("url_export/PMID36635386_studies_export.tsv")
export_data <- export_data[grepl(" European", export_data$discoverySampleAncestry),]
urls = paste0(export_data$summaryStatistics, "/harmonised/", export_data$accessionId, ".h.tsv.gz")
write.table(urls, "url_txt/1400met_urls.txt", row.names = F, col.names = F, quote = F)

#https://www.ebi.ac.uk/gwas/publications/38448586
export_data <- data.table::fread("url_export/PMID38448586_studies_export.tsv")
urls = paste0(export_data$summaryStatistics, "/harmonised/", export_data$accessionId, ".h.tsv.gz")
write.table(urls, "url_txt/233met_urls.txt", row.names = F, col.names = F, quote = F)

#https://www.ebi.ac.uk/gwas/publications/33437055
export_data <- data.table::fread("url_export/PMID33437055_studies_export.tsv")
urls = paste0(export_data$summaryStatistics, "/harmonised/33437055-", export_data$accessionId, "-EFO_0004725.h.tsv.gz")
write.table(urls, "url_txt/338met_urls.txt", row.names = F, col.names = F, quote = F)
