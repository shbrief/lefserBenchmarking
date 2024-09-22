#' Create input for Venn Diagram plot by eulerr package
#' 
#' @importFrom eulerr euler
#'
#' @param lefse_taxa A character vector.
#' @param lefser_taxa A character vector.
#' @param mm_taxa A character vector.
#'
#' @return An `euler` object ready for `plot` function.
#' 
#' @export
threeVennDiagram <- function(lefse_taxa, lefser_taxa, mm_taxa) {
    ABC <- length(intersect(lefse_taxa, lefser_taxa) %>% intersect(., mm_taxa))
    AB <- length(intersect(lefse_taxa, lefser_taxa)) - ABC
    AC <- length(intersect(lefse_taxa, mm_taxa)) - ABC
    BC <- length(intersect(lefser_taxa, mm_taxa)) - ABC
    
    fit <- eulerr::euler(c(
        "LEfSe" = length(lefse_taxa) -AB - AC -ABC, # A
        "lefser" = length(lefser_taxa) -AB -BC -ABC, # B
        "MM" = length(mm_taxa) -AC -BC -ABC, # C
        "LEfSe&lefser" = AB,
        "LEfSe&MM" = AC,
        "lefser&MM" = BC,
        "LEfSe&lefser&MM" = ABC
    ))
    
    return(fit)
}



#' Create input for Venn Diagram plot by eulerr package
#' 
#' @importFrom eulerr euler
#'
#' @param input1 A character vector.
#' @param input2 A character vector.
#' @param name1 A character(1).
#' @param name2 A character(1).
#'
#' @return An `euler` object ready for `plot` function.
#' 
#' @export
twoVennDiagram <- function(input1, input2, name1 = "A", name2 = "B") {
    
    AB <- length(intersect(input1, input2))
    
    summary <- c(length(input1) -AB, length(input2) -AB, AB)
    names(summary) <- c(name1, name2, paste0(name1, "&", name2))
    
    fit <- eulerr::euler(summary)
    
    return(fit)
}