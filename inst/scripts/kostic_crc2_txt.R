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
    select(DIAGNOSIS, X.SampleID) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('helper_col') %>%
    mutate(across(.cols = everything(), .fns = ~as.character(.x))) %>%
    set_colnames(paste0('col1', seq_along(.)))

## Combine counts and sample_metadata in a single data frame
data <- bind_rows(sm, counts)
colnames(data) <- NULL

## Remove arbitrary `OTU__` part
data[1] <- lapply(data[1], function(x) {
    gsub("OTU__[0-9]+\\|", "", x)}) %>% unlist

## Export to txt file
write.table(
    data, file = file.path(projDir, "kostic_crc2.txt"),
    sep = '\t', row.names = FALSE, col.names = FALSE,
    quote = FALSE
)
