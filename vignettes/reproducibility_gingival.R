suppressPackageStartupMessages({
    library(microbiomeMarker)
    library(MicrobiomeBenchmarkData)
    library(MicrobiomeBenchmarkDataAnalyses)
    library(lefser)
    library(dplyr)
    library(purrr)
    library(stringr)
    library(phyloseq)
    library(mia)
})

## Parameters
kwth <- 0.01
wth <- 0.01
ldath <- 3


##### LEfSe ----------------------------

## LEfSe results from Conda
LEfSeDatDir <- "~/Projects/MicrobiomeBenchmarkDataLefse/inst/extdata"
col_names <- c('feature', 'log_hi_class_avg', 'class', 'lefse_conda_LDA', 'pval')
fnames <- list.files(LEfSeDatDir) %>% .[grep(".res", .)]
lda3_fnames <- fnames[grep("lda3", fnames)]

## A list of 10 LEfSe iterations
LEfSeResAll <- as.list(vector(length = 10))
for (i in seq_along(lda3_fnames)) {
    res <- readr::read_tsv(
        file.path(LEfSeDatDir, lda3_fnames[i]), col_names = FALSE, show_col_types = FALSE) %>% 
        magrittr::set_colnames(col_names) %>% 
        filter(!is.na(lefse_conda_LDA)) %>% 
        mutate(app_name = 'lefse_conda')
    
    res$feature <- lapply(res$feature, function(x) {
        strsplit(x, "\\.") %>% unlist %>% tail(., 1) %>% 
            gsub("t__|c__|f__|s__|p__|o__|g__", "", .)}) %>% unlist
    
    LEfSeResAll[[i]] <- res
}

## Combine LDA scores 
LEfSe_iter_tb <- LEfSeResAll[[1]][c("feature", "lefse_conda_LDA")] %>% 
    magrittr::set_colnames(c("feature", paste0("LDA", 1)))

for (i in 2:10) {
    sub <- LEfSeResAll[[i]][c("feature", "lefse_conda_LDA")] %>% 
        magrittr::set_colnames(c("feature", paste0("LDA", i)))
    LEfSe_iter_tb <- dplyr::full_join(x = LEfSe_iter_tb, y = sub, by = "feature")
}

LEfSe_iter_tb_nonNA <- LEfSe_iter_tb
LEfSe_iter_tb_nonNA[is.na(LEfSe_iter_tb_nonNA)] <- 0

## Calculate mean and sd of iterations
LEfSe_iter_tb_nonNA$mean <- apply(LEfSe_iter_tb_nonNA[2:11], 1, mean)
LEfSe_iter_tb_nonNA$sd <- apply(LEfSe_iter_tb_nonNA[2:11], 1, sd)

## Save
write.csv(LEfSe_iter_tb_nonNA, 
          "~/Projects/lefserBenchmarking/data/gingival_Outputs_Iterations/LEfSe_10iters.csv",
          row.names = FALSE)



##### lefser --------------------------
## Data
dat_name <- 'HMP_2012_16S_gingival_V35'
conditions_col <- 'body_subsite'
conditions <- c(condB = 'subgingival_plaque', condA = 'supragingival_plaque')
tse <- getBenchmarkData(dat_name, dryrun = FALSE)[[1]]
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
rankNames <- colnames(rowData(tse_subset))
rankNames <- stringr::str_replace(rankNames, "superkingdom", "kingdom")
colnames(rowData(tse_subset)) <- rankNames

lefserInput <- relativeAb(tse_subset)
colData(lefserInput)$body_subsite <- factor(
    colData(lefserInput)$body_subsite,
    levels = c("supragingival_plaque", "subgingival_plaque")
)

## A list of 10 lefser iterations
lefserResAll <- as.list(vector(length = 10))

for (i in 1:10) {
    set.seed(i)
    res <- lefser(
        relab = lefserInput,
        kruskal.threshold = kwth,
        wilcox.threshold = wth,
        lda.threshold = ldath,
        groupCol = "body_subsite",
        blockCol = NULL
    )
    lefserResAll[[i]] <- res
}

## Combine LDA scores
lefser_iter_tb <- lefserResAll[[1]][c("features", "scores")] %>% 
    magrittr::set_colnames(c("feature", paste0("LDA", 1)))

for (i in 2:10) {
    sub <- lefserResAll[[i]][c("features", "scores")] %>% 
        magrittr::set_colnames(c("feature", paste0("LDA", i)))
    lefser_iter_tb <- dplyr::full_join(x = lefser_iter_tb, y = sub, by = "feature")
}

lefser_iter_tb_nonNA <- lefser_iter_tb
lefser_iter_tb_nonNA[is.na(lefser_iter_tb_nonNA)] <- 0

## Calculate mean and sd of iterations
lefser_iter_tb_nonNA$mean <- apply(lefser_iter_tb_nonNA[2:11], 1, mean)
lefser_iter_tb_nonNA$sd <- apply(lefser_iter_tb_nonNA[2:11], 1, sd)

## Save
write.csv(lefser_iter_tb_nonNA, 
          "~/Projects/lefserBenchmarking/data/gingival_Outputs_Iterations/lefser_10iters.csv",
          row.names = FALSE)


##### microbiomeMarker ---------------------------
## Data
ps <- convertToPhyloseq(tse_subset)
colnames(tax_table(ps)) <- str_to_sentence(colnames(tax_table(ps)))

## A list of 10 LEfSe iterations
mmResAll <- as.list(vector(length = 10))
for (i in 1:10) {
    set.seed(i)
    res <- run_lefse(
        ps = ps,
        group = "body_subsite",
        taxa_rank = "none",
        norm = "CPM", # This is the relative abundance but with a total ofset to 1e6
        kw_cutoff = kwth,
        wilcoxon_cutoff = wth,
        lda_cutoff = ldath,
        multigrp_strat = FALSE,
        bootstrap_n = (30 + i*10)
    )
    
    mmResTb <- marker_table(res)
    class(mmResTb) <- NULL
    mmResTb <- as_tibble(mmResTb)[c("feature", "ef_lda")]
    
    mmResAll[[i]] <- mmResTb
}

## Combine LDA scores 
mm_iter_tb <- mmResAll[[1]][c("feature", "ef_lda")] %>% 
    magrittr::set_colnames(c("feature", paste0("LDA", 1)))

for (i in 2:10) {
    sub <- mmResAll[[i]][c("feature", "ef_lda")] %>% 
        magrittr::set_colnames(c("feature", paste0("LDA", i)))
    mm_iter_tb <- dplyr::full_join(x = mm_iter_tb, y = sub, by = "feature")
}

mm_iter_tb_nonNA <- mm_iter_tb
mm_iter_tb_nonNA[is.na(mm_iter_tb_nonNA)] <- 0

## Calculate mean and sd of iterations
mm_iter_tb_nonNA$mean <- apply(mm_iter_tb_nonNA[2:11], 1, mean)
mm_iter_tb_nonNA$sd <- apply(mm_iter_tb_nonNA[2:11], 1, sd)

## Save
write.csv(mm_iter_tb_nonNA, 
          "~/Projects/lefserBenchmarking/data/gingival_Outputs_Iterations/microbiomeMarker_10iters.csv",
          row.names = FALSE)


