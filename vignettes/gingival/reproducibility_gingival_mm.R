suppressPackageStartupMessages({
    library(microbiomeMarker)
    library(lefser)
    library(dplyr)
    library(purrr)
    library(stringr)
    library(phyloseq)
    library(mia)
})

##### Prepare input data -------------------------------------
tse_subset <- readRDS("data/gingival_dataset_for_benchmarking.rds")
rankNames <- colnames(rowData(tse_subset))
rankNames <- stringr::str_replace(rankNames, "superkingdom", "kingdom")
colnames(rowData(tse_subset)) <- rankNames
colData(tse_subset)$body_subsite <- factor(
    colData(tse_subset)$body_subsite,
    levels = c("supragingival_plaque", "subgingival_plaque")
)

tse_subset <- relativeAb(tse_subset)
ps <- convertToPhyloseq(tse_subset)
colnames(tax_table(ps)) <- str_to_sentence(colnames(tax_table(ps)))


##### microbiomeMarker: single run -------------------------------------------
## Parameters
kwth <- 0.05
wth <- 0.05
ldath <- 2
# kwth <- 0.01
# wth <- 0.01
# ldath <- 3

set.seed(1982)
res <- run_lefse(
    ps = ps,
    group = "body_subsite",
    taxa_rank = "none",
    norm = "CPM", # This is the relative abundance but with a total of set to 1e6
    kw_cutoff = kwth,
    wilcoxon_cutoff = wth,
    lda_cutoff = ldath,
    multigrp_strat = FALSE
)

fnames <- paste0("gingival_a", kwth, "_w", wth, "_lda", ldath, "_mm.rds")
saveRDS(res, file.path("data", fnames))


# ##### 10 iterations of microbiomeMarker ------------------------------
# ## Parameters
# kwth <- 0.01
# wth <- 0.01
# ldath <- 3
# 
# mmResAll <- as.list(vector(length = 10))
# for (i in 1:10) {
#     set.seed(i)
#     res <- run_lefse(
#         ps = ps,
#         group = "body_subsite",
#         taxa_rank = "none",
#         norm = "CPM", # This is the relative abundance but with a total ofset to 1e6
#         kw_cutoff = kwth,
#         wilcoxon_cutoff = wth,
#         lda_cutoff = ldath,
#         multigrp_strat = FALSE,
#         bootstrap_n = (30 + i*10)
#     )
# 
#     mmResTb <- marker_table(res)
#     class(mmResTb) <- NULL
#     mmResTb <- as_tibble(mmResTb)[c("feature", "ef_lda")]
# 
#     mmResAll[[i]] <- mmResTb
# }
# 
# ## Combine LDA scores
# mm_iter_tb <- mmResAll[[1]][c("feature", "ef_lda")] %>%
#     magrittr::set_colnames(c("feature", paste0("LDA", 1)))
# 
# for (i in 2:10) {
#     sub <- mmResAll[[i]][c("feature", "ef_lda")] %>%
#         magrittr::set_colnames(c("feature", paste0("LDA", i)))
#     mm_iter_tb <- dplyr::full_join(x = mm_iter_tb, y = sub, by = "feature")
# }
# 
# mm_iter_tb_nonNA <- mm_iter_tb
# mm_iter_tb_nonNA[is.na(mm_iter_tb_nonNA)] <- 0
# 
# ## Calculate mean and sd of iterations
# mm_iter_tb_nonNA$mean <- apply(mm_iter_tb_nonNA[2:11], 1, mean)
# mm_iter_tb_nonNA$sd <- apply(mm_iter_tb_nonNA[2:11], 1, sd)
# 
# ## Save
# write.csv(mm_iter_tb_nonNA,
#           "~/Projects/lefserBenchmarking/data/gingival_Outputs_Iterations/microbiomeMarker_10iters.csv",
#           row.names = FALSE)
