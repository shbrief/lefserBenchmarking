.empty_func <- function(x) {
    if(length(x) == 0) {
        return(NA)
    } else {return(x)}
}

#' Split the taxonomic label in character into the tax_table format
#'
#' @param x A character vector of taxonomic names
#' @param prefix Logical. If it's `FALSE`, there is a prefix showing the
#' taxonomic levels.
#'
#' @export
taxNameToTable <- function(x, prefix = TRUE) {

    classes <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    taxTb <- matrix(nrow = length(x), ncol = length(classes))
    colnames(taxTb) <- classes
    rownames(taxTb) <- x

    for (i in seq_along(x)) {
        if (isTRUE(prefix)) {
            y <- strsplit(x[i], "\\|") %>% unlist

            taxTb[i, "Kingdom"] <- y[grep("k__", y)] %>% gsub("k__", "", .) %>% .empty_func()
            taxTb[i, "Phylum"] <- y[grep("p__", y)] %>% gsub("p__", "", .) %>% .empty_func()
            taxTb[i, "Class"] <- y[grep("c__", y)] %>% gsub("c__", "", .) %>% .empty_func()
            taxTb[i, "Order"] <- y[grep("o__", y)] %>% gsub("o__", "", .) %>% .empty_func()
            taxTb[i, "Family"] <- y[grep("f__", y)] %>% gsub("f__", "", .) %>% .empty_func()
            taxTb[i, "Genus"] <- y[grep("g__", y)] %>% gsub("g__", "", .) %>% .empty_func()
            taxTb[i, "Species"] <- y[grep("s__", y)] %>% gsub("s__", "", .) %>% .empty_func()
        } else {
            z <- strsplit(x[i], "\\|") %>% unlist
            l <- length(z)
            z_filled <- c(z, rep(NA, abs(length(classes) - l)))
            taxTb[i,] <- z_filled
        }
    }
    return(taxTb)
}
