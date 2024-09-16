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

lefserInput <- relativeAb(tse_subset)


##### lefser -------------------------------------------
## Parameters
kwth <- 0.05
wth <- 0.05
ldath <- 2
# kwth <- 0.01
# wth <- 0.01
# ldath <- 3

set.seed(1982)
res <- lefser(
    relab = lefserInput,
    kruskal.threshold = kwth,
    wilcox.threshold = wth,
    lda.threshold = ldath,
    groupCol = "body_subsite",
    blockCol = NULL
)

fnames <- paste0("gingival_a", kwth, "_w", wth, "_lda", ldath, "_lefser.rds")
saveRDS(res, file.path("data", fnames))


##### lefser: terminal nodes -------------------------------------------
## Parameters
kwth <- 0.01
wth <- 0.01
ldath <- 3

## Only with the terminal nodes
terminal_nodes <- get_terminal_nodes(rownames(tse_subset))
lefserInput2 <- tse_subset[terminal_nodes, ] %>% relativeAb(.)

set.seed(1982)
res2 <- lefser(
    relab = lefserInput2,
    kruskal.threshold = kwth,
    wilcox.threshold = wth,
    lda.threshold = ldath,
    groupCol = "body_subsite",
    blockCol = NULL
)
saveRDS(res2, "data/lefser_terminal_nodes.rds")

