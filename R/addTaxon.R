#' Add prior information using 'features' as key
#' 
#' @param x A data frame with the key column labeled as `features`
#' @param prior_info 
#' @param reformat Logical. Under the default (`TRUE`), any arbitrary 
#' delimiters will be removed and the rest will be concatenated. 
#'
#' @export
addTaxon <- function(x, prior_info, reformat = TRUE) {
    
    ## Set the feature level based on average LDA scores
    if ("mean" %in% colnames(x)) {
        lvs <- x$features[order(x$mean, decreasing = TRUE)]
        x$features <- factor(x$features, level = lvs)
    }
    
    if (isTRUE(reformat)) {
        x$features <- x$features %>% gsub("_", "", .) %>% gsub("\\.", "", .)
    }
    
    x_taxon <- x %>% 
        # dplyr::rename(features = feature) %>%
        dplyr::left_join(., prior_info, by = "features")
    
    return(x_taxon)
}