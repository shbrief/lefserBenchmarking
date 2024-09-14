#' Function to add taxonomic level prefixes
#' 
#' @param taxonomic_names A character vector
#' @param current_delim A character(1). Regular expression of the current delimiter
#' 
#' @examples
#' taxonomic_name <- "Bacteria.Firmicutes.Bacilli.Lactobacillales.Streptococcaceae"
#' addTaxonomicPrefixes(taxonomic_name, "\\.")
#' 
#' @export
addTaxonomicPrefixes <- function(taxonomic_names, current_delim) {
    # Define the taxonomic levels
    levels <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
    
    res <- vector(mode = "character", length = length(taxonomic_names))
    for (i in seq_along(taxonomic_names)) {
        # Split the taxonomic name into its components
        components <- strsplit(taxonomic_names[i], current_delim)[[1]]
        
        # Ensure we don't exceed the number of defined levels
        num_components <- min(length(components), length(levels))
        
        # Add prefixes to each component
        prefixed_components <- paste0(levels[1:num_components], components[1:num_components])
        
        # Join the components back together using '|' as separator
        result <- paste(prefixed_components, collapse = "|")
        
        res[i] <- result
    }
    
    return(res)
}



#' 
#' @export
getMostSpecifiTaxa <- function(taxonomic_names, current_delim) {
    res <- taxonomic_names %>%
        addTaxonomicPrefixes(., current_delim = current_delim) %>%
        strsplit(., "\\|") %>% 
        lapply(., function(x) {tail(x, 1)}) %>% unlist
    return(res)
}
