tse_subset <- readRDS("data/gingival_dataset_for_benchmarking.rds")
tse_subset <- relativeAb(tse_subset) # relative abundance

## Count part of the table
counts <- tse_subset |>
    assay() |>
    as.data.frame() %>%
    tibble::rownames_to_column('feature') |>
    mutate(across(.cols = everything(), .fns = ~as.character(.x))) |>
    mutate(feature = gsub("\\.", "_", feature)) |>
    mutate(feature = gsub("_", "", feature))
colnames(counts) <- paste0("col", seq_along(counts))

## Metadata part of the table
sm <- tse_subset |>
    colData() |>
    as.data.frame() |>
    tibble::rownames_to_column('Sample') |>
    # select(body_subsite, gender, Sample) |>
    select(body_subsite, Sample) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column('helper_col') |>
    mutate(across(.cols = everything(), .fns = ~ as.character(.x)))
colnames(sm) <- paste0('col', seq_along(sm))

dat <- bind_rows(sm, counts)
colnames(dat) <- NULL

## Save
write.table(
    dat, 'data/gingivalplaque_ra.txt', sep = '\t', row.names = FALSE, col.names = FALSE,
    quote = FALSE
)
