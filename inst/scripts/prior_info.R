library(MicrobiomeBenchmarkData)

dat_name <- 'HMP_2012_16S_gingival_V35'
tse <- getBenchmarkData(dat_name, dryrun = FALSE)[[1]]
row_data <- as.data.frame(rowData(tse))
prior_info <- row_data[, c('genus', 'taxon_annotation')]
prior_info$taxon_name <- rownames(row_data)
prior_info$new_names <- paste0(prior_info$taxon_name, '|', prior_info$genus)
prior_info$features <- gsub("_", "", prior_info$taxon_name) %>% gsub("\\.", "", .)
prior_info <- 
    dplyr::relocate(prior_info, taxon_name, features, new_names, genus, taxon_annotation)
write.csv(prior_info, "~/Projects/lefserBenchmarking/inst/extdata/prior_info.csv", row.names = FALSE)