#'
#' @param taxTable A `taxonomyTable` object from phyloseq.
#'
#' @examples
#' data(kostic_crc_small)
#' taxTableToName(tax_table(kostic_crc_small))

taxTableToName <- function(taxTable) {
    
    if (class(taxTable) != "taxonomyTable") {
        msg <- "The input should be taxonomyTable object from phyloseq."
        stop(msg)
    }
    
    dat <- as.data.frame(taxTable)
    ori_colnames <- colnames(dat)
    colnames(dat) <- colnames(dat) %>% 
        substr(., 1, 1) %>% 
        tolower()
    dat <- tibble::rownames_to_column(dat, "OTU")

    taxnamed <- OmicsMLRepoR::getNarrowMetaTb(
        dat,
        newCol = "taxname",
        targetCols = colnames(dat),
        sep = "__",
        delim = "|",
        remove = TRUE,
        na.rm = TRUE,
        sort = FALSE
    ) 
    
    colnames(dat) <- ori_colnames # put original column names back
    res <- cbind(taxnamed, dat) %>% tibble::column_to_rownames(var = "taxname")
    return(res)
}