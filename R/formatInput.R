#' Format metagenomic data for different biomarker discovery tools
#'
#' Input format conversion across `SummarizedExperiment` and `phyloseq` objects.
#'
#' @import phyloseq
#' @import SummarizedExperiment
#' @importFrom mia makePhyloseqFromTreeSE
#'
#' @param dat A microbiome abundance data. It can be a data frame (an input
#' format for LEfSe's `run_lefse.py` function), `SummarizedExperiment`, or
#' `phyloseq` object.
#' @param format_to The target format. Available options are
#' `c("SummarizedExperiment", "phyloseq")`.
#'
#' @examples
#' formatInput(kostic_crc_small, format_to = "SummarizedExperiment")
#'
#' @export
formatInput <- function(dat, format_to) {

    ##### phyloseq to SE ----------------
    if (class(dat) == "phyloseq") {

        ##### Deconstruct
        assay <- unclass(otu_table(dat))
        coldata <- as(sample_data(dat), "data.frame")
        rowdata <- tax_table(dat) %>% taxTableToName()
        rownames(assay) <- rownames(rowdata) # UPDATE to link through key column

        ##### Reconstruct
        se <- SummarizedExperiment(
            assays = list(exprs = assay), colData = coldata, rowData = rowdata
        )
        return(se)
    }

    ##### SE to phyloseq -----------------
    if (class(dat) == "SummarizedExperiment") {

        ##### Deconstruct
        assay <- assays(dat)[[1]] # UPDATE to default the first assay (not hard-code)
        coldata <- colData(dat)
        rowdata <- rowData(dat)
        if (ncol(rowdata) == 0) {
            rowdata <- taxNameToTable(rownames(rowdata))
        }

        ##### Reconstruct
        ## Arbitrary OTU
        rownames(assay) <- seq_len(nrow(assay))
        rownames(rowdata) <- seq_len(nrow(rowdata))

        # ## Format rownames <<<<<<<<< Not sure the format requirement for otu_table
        # all_prefix <- "k__|p__|c__|o__|f__|g__|s__"
        # rownames(assay) <- gsub(all_prefix, "", rownames(assay))
        # rownames(rowdata) <- gsub(all_prefix, "", rownames(rowdata))

        ## Assemble into SummarizedExperiment
        se <- SummarizedExperiment(
            assays = list(exprs = assay), colData = coldata, rowData = rowdata
        )

        ## Use the mia package for conversion
        ps <- mia::makePhyloseqFromTreeSE(se, assay.type = "exprs")
        return(ps)
    }

}
