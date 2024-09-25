#' Get the LEfSe iteration summary table 
#'
#' @return A data frame with '`iter`+3' columns. `iter` columns are results from
#' each bootstrap iterations and three columns are `feature`, `mean`, and `sd`. 
#'
#' @export
getLEfSeIterSummary <- function(LEfSeDatDir, iter = 10, reference = NA) {
    
    ## LEfSe results from Docker
    col_names <- c('feature', 'log_hi_class_avg', 'class', 'LDA', 'pval')
    fnames <- list.files(LEfSeDatDir) %>% .[grep(".res", .)]
    
    ## Order by the number of iterations
    bp_values <- as.numeric(sub(".*bp(\\d+)\\.res", "\\1", fnames))
    ordered_fnames <- fnames[order(bp_values)]
    
    ## A list of `iter` LEfSe iterations (default = 10)
    LEfSeResAll <- as.list(vector(length = iter))
    for (i in seq_along(ordered_fnames)) {
        
        iter <- sub(".*bp(\\d+)\\.res", "\\1", ordered_fnames[i])
        res <- readr::read_tsv(
            file.path(LEfSeDatDir, ordered_fnames[i]), col_names = FALSE, show_col_types = FALSE) %>% 
            magrittr::set_colnames(col_names) %>% 
            filter(!is.na(LDA))
        
        ## Mark the enrichment direction
        if (is.na(reference)) {
            lvs <- paste(unique(res$class), collapse = ", ")
            msg <- paste("`reference` level is missing. Specify from followings:\n",
                         lvs)
            stop(msg)
        }
        res$LDA <- ifelse(res$class == reference, res$LDA*(-1), res$LDA)
        
        # res$feature <- lapply(res$feature, function(x) {
        #     strsplit(x, "\\.") %>% unlist %>% tail(., 1) %>% 
        #         gsub("k__|p__|c__|o__|f__|g__|s__", "", .)}) %>% unlist
        
        LEfSeResAll[[i]] <- res
        names(LEfSeResAll)[i] <- paste0("iter_", iter)
    }
    
    ## Combine LDA scores 
    LEfSe_iter_tb <- LEfSeResAll[[1]][c("feature", "LDA")] %>% 
        magrittr::set_colnames(c("feature", paste0("LDA_", names(LEfSeResAll)[1])))
    
    for (i in 2:length(LEfSeResAll)) {
        sub <- LEfSeResAll[[i]][c("feature", "LDA")] %>% 
            magrittr::set_colnames(c("feature", paste0("LDA_", names(LEfSeResAll)[i])))
        LEfSe_iter_tb <- dplyr::full_join(x = LEfSe_iter_tb, y = sub, by = "feature")
    }
    
    LEfSe_iter_tb_nonNA <- LEfSe_iter_tb
    LEfSe_iter_tb_nonNA[is.na(LEfSe_iter_tb_nonNA)] <- 0
    
    ## Calculate mean and sd of iterations
    LEfSe_iter_tb_nonNA$mean <- apply(LEfSe_iter_tb_nonNA[-1], 1, mean)
    LEfSe_iter_tb_nonNA$sd <- apply(LEfSe_iter_tb_nonNA[-1], 1, sd)
    
    return(LEfSe_iter_tb_nonNA)
}