library(MicrobiomeBenchmarkData)

dat_name <- 'HMP_2012_16S_gingival_V35'
tse <- getBenchmarkData(dat_name, dryrun = FALSE)[[1]]
row_data <- as.data.frame(rowData(tse))

## Function to get the most specific non-NA taxonomic name
get_most_specific_taxa <- function(row) {
    taxa_levels <- c("genus", "family", "order", "class", "phylum", "superkingdom")
    for (level in taxa_levels) {
        if (!is.na(row[level]) && row[level] != "<NA>") {
            return(paste(substring(level, 1, 1), row[level], sep = "__"))
        }
    }
    return("No valid taxonomic information")
}
row_data$most_specific_taxa <- apply(row_data, 1, get_most_specific_taxa)

prior_info <- row_data[, c("most_specific_taxa", "taxon_annotation")]
prior_info$taxon_name <- rownames(row_data)
prior_info$new_names <- paste0(prior_info$taxon_name, '|', prior_info$most_specific_taxa)
prior_info$features <- gsub("_", "", prior_info$taxon_name) %>% gsub("\\.", "", .)
prior_info <- 
    dplyr::relocate(prior_info, taxon_name, features, new_names, most_specific_taxa, taxon_annotation)
write.csv(prior_info, "~/Projects/lefserBenchmarking/inst/extdata/prior_info.csv", row.names = FALSE)
