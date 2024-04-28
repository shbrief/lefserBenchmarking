library(microbiomeMarker)
library(lefserBenchmarking)

## Load data
data(kostic_crc)
kostic_crc_small <- phyloseq::subset_taxa(
    kostic_crc,
    Phylum %in% c("Firmicutes")
)

## Convert phyloseq object to SummarizedExperiment object
kostic_se <- formatInput(kostic_crc_small, format_to = "SummarizedExperiment")
kostic_se$DIAGNOSIS <- as.factor(kostic_se$DIAGNOSIS)

## Convert SummarizedExperiment to a LEfSe input table
exprs <- assay(kostic_se, i = 1L)
datas <- data.frame(
    condition = as.factor(kostic_se$DIAGNOSIS),
    subjectID = kostic_se$X.SampleID,
    t(exprs)
)

text <- mapply(function(x, y) {paste0(x, "\t", y)}, 
               x = names(datas), 
               y = sapply(datas, paste0, collapse = "\t"))

inputDir <- "~/Projects/lefserBenchmarking/LEfSe_Inputs"
writeLines(text, 
           con = file(file.path(inputDir, "kostic_crc.txt")))