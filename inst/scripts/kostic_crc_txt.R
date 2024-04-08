## Modifed https://github.com/sdgamboa/lefse_comparison/blob/main/get_dataset.R
## Prepare LEfSe input table

suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(microbiomeMarker)
    library(dplyr)
    library(tibble)
    library(magrittr)
    library(lefserBenchmark)
})

projDir <- "~/Projects/lefserBenchmark/data/LEfSe_Inputs"

## Covert phylose object to SummarizedExperiment object
data(kostic_crc)
kostic_crc_small <- phyloseq::subset_taxa(
    kostic_crc,
    Phylum %in% c("Firmicutes")
)

se <- formatInput(kostic_crc_small, format_to = "SummarizedExperiment")
se$DIAGNOSIS <- as.factor(se$DIAGNOSIS)

## Example datset from lefser package as SummarizedExperiment
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
    data, file = file.path(projDir, "kostic_crc.txt"),
    sep = '\t', row.names = FALSE, col.names = FALSE,
    quote = FALSE
)
