library(MicrobiomeBenchmarkData)
library(MicrobiomeBenchmarkDataAnalyses)
library(dplyr)
library(purrr)
library(tidySummarizedExperiment)

## Gingival data from the MicrobiomeBenchmarkData package
dat_name <- 'HMP_2012_16S_gingival_V35'
conditions_col <- 'body_subsite'
conditions <- c(condB = 'subgingival_plaque', condA = 'supragingival_plaque')
tse <- getBenchmarkData(dat_name, dryrun = FALSE)[[1]]

## Extract information
col_data <- tse |>
    colData() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_name") |>
    as_tibble()
subjects <- col_data |>
    pull(subject_id) |>
    unique()
sample_names <- vector("list", length(subjects))
names(sample_names) <- subjects

## Filter only the paired samples from a single visit
for (i in seq_along(subjects))  {
    current_subject <- subjects[i]
    sub_dat <- col_data |>
        filter(subject_id == current_subject) |>
        slice_max(order_by = visit_number, with_ties = TRUE, n = 1)
    if (nrow(sub_dat) < 2) {
        next
    }
    lgl_vct <- all(sort(sub_dat[[conditions_col]]) == conditions)
    if (isFALSE(lgl_vct)) {
        next
    }
    sample_names[[i]] <- sub_dat
}
sample_names <- discard(sample_names, is.null)
col_data_subset <- bind_rows(sample_names)
selected_samples <- col_data_subset |>
    pull(sample_name)
tse_subset <- tse[, selected_samples]
tse_subset <- filterTaxa(tse_subset)

## Save
saveRDS(tse_subset, "data/gingival_dataset_for_benchmarking.rds")
