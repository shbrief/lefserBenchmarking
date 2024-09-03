## From https://github.com/sdgamboa/lefse_comparison/blob/main/get_dataset.R
## Convert SummarizedExperiment object to LEfSe input table

suppressPackageStartupMessages({
    library(lefser)
    library(SummarizedExperiment)
    library(dplyr)
    library(tibble)
    library(magrittr)
})

projDir <- "~/Projects/lefserBenchmark/data/LEfSe_Inputs"

## Example dataset from the lefser package as SummarizedExperiment
data(zeller14)
se <- zeller14[, zeller14$study_condition != "adenoma"]
rownames(se) <- sub('^.+([a-z]__.+$)', '\\1', rownames(se))

## Matrix with counts
counts <- assay(se) %>%
    as.data.frame() %>%
    rownames_to_column('features') %>%
    mutate(across(.cols = everything(), .fns = ~as.character(.x))) %>%
    set_colnames(paste0('col1', seq_along(.)))

## Sample metadata
sm <- colData(se) %>%
    as.data.frame() %>%
    rownames_to_column('Sample') %>%
    select(study_condition, age_category, Sample) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('helper_col') %>%
    mutate(across(.cols = everything(), .fns = ~as.character(.x))) %>%
    set_colnames(paste0('col1', seq_along(.)))

## Combine counts and sample_metadata in a single data frame
data <- bind_rows(sm, counts)
colnames(data) <- NULL

## Export to txt file
write.table(
    data, file.path(projDir, "zeller14.txt"),
    sep = '\t', row.names = FALSE, col.names = FALSE,
    quote = FALSE
)
