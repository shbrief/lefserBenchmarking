# Mice data --------------------------
suppressPackageStartupMessages({
    library(lefser)
    library(microbiomeMarker)
    library(dplyr)
    library(ggplot2)
    library(lefserBenchmarking)
    library(grid)
    library(VennDiagram)
})

projDir <- "~/Projects/lefserBenchmarking"
mice <- readr::read_table(file.path(projDir, "data/LEfSe_Inputs/mice_rag2.txt"), 
                          col_names = FALSE)

## Put together as SummarizedExperiment
relative_ab <- mice[-1,] %>% tibble::column_to_rownames("X1") 
relative_ab <- apply(relative_ab, 2, as.numeric)
rownames(relative_ab) <- mice[["X1"]][-1]
colData <- as(mice[1, -1], "data.frame") %>% t 
colnames(colData) <- "Genotype"
rowData <- taxNameToTable(rownames(relative_ab), FALSE)
mice_se <- SummarizedExperiment(
    assays = list(relative_abundance = relative_ab), 
    colData = colData, 
    rowData = rowData
)

## Run lefser -----------
set.seed(1982)
res <- lefser(relativeAb(mice_se), 
              groupCol = "Genotype", 
              kruskal.threshold = 0.01,
              lda.threshold = 2)

## Clean the lefser result table
lefser_output <- res %>% 
    mutate(app_name = 'lefser') %>% 
    arrange(scores) %>% 
    dplyr::rename(lefser_LDA = scores)
lefser_output$feature <- lapply(lefser_output$feature, function(x) {
    strsplit(x, "\\|") %>% unlist %>% tail(., 1)}) %>% unlist

## Run microbiomeMarker (MM) ---------
mice_ps <- mia::convertToPhyloseq(mice_se, 
                                  assay.type = "relative_abundance")
set.seed(1982)
mm_lefse <- run_lefse(
    mice_ps,
    wilcoxon_cutoff = 0.01,
    group = "Genotype",
    multigrp_strat = TRUE,
    lda_cutoff = 2
)

## Clean the MM result table
mm_output <- data.frame(
    feature = marker_table(mm_lefse)$feature,
    mm_LDA = marker_table(mm_lefse)$ef_lda,
    app_name = "microbiomeMarker") %>% arrange(mm_LDA)
mm_output$feature <- lapply(mm_output$feature, function(x) {
    strsplit(x, "\\|") %>% unlist %>% tail(., 1)}) %>% unlist

## LEfSe result from Docker --------------
LEfSeZellerDir <- "~/Projects/lefserBenchmarking/data/LEfSe_Outputs/mice_rag2"
col_names <- c('feature', 'log_hi_class_avg', 'class', 'lefse_docker_LDA', 'pval')
lefse_docker <- readr::read_tsv(
    file.path(LEfSeZellerDir, "mice_rag2.res"),
    col_names = FALSE, show_col_types = FALSE) %>% 
    magrittr::set_colnames(col_names) %>% 
    filter(!is.na(lefse_docker_LDA)) %>% 
    mutate(lefse_galaxy_LDA = ifelse(class == "control", -lefse_docker_LDA, lefse_docker_LDA),
           app_name = 'lefse_docker')

lefse_docker$feature <- lapply(lefse_docker$feature, function(x) {
    strsplit(x, "\\.") %>% unlist %>% tail(., 1) %>% gsub("t__|c__|f__|s__|p__|o__|g__", "", .)}) %>% unlist
