#' @param x A data frame with the key column labeled as `features`
#'
#' @export
addTaxon <- function(x, prior_info) {
    
    ## Set the feature level based on average LDA scores
    if ("mean" %in% colnames(x)) {
        lvs <- x$features[order(x$mean, decreasing = TRUE)]
        x$features <- factor(x$features, level = lvs)
    }
    
    x$features <- x$features %>% gsub("_", "", .) %>% gsub("\\.", "", .)
    x_taxon <- x %>% 
        # dplyr::rename(features = feature) %>%
        dplyr::left_join(., prior_info, by = "features")
    
    return(x_taxon)
}