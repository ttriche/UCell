#' Calculate UCell scores on SingleCellExperiment objects
#'
#' A wrapper for the function \code{\link{ScoreSignatures_UCell}} that interacts directly with SingleCellExperiment (sce) objects. (This could be an S4 generic)
#'
#' @param sce       A SingleCellExperiment object
#' @param assay     The sce object assay where the data is to be found
#' @param features  A list of signatures, for example: \code{list( Tcell_signature = c("CD2","CD3E","CD3D"), Myeloid_signature = c("SPI1","FCER1G","CSF1R"))}
#'     You can also specify positive and negative gene sets by adding a + or - sign to genes in the signature; see an example below
#' @param ... Additional parameters to be passed on to \code{\link{ScoreSignatures_UCell}}
#' @return the SingleCellExperiment, with UCell scores added to altExp field
#' @examples
#' ## Not run:
#' library(UCell)
#' my.matrix <- UCell::sample.matrix
#' my.sce <- SingleCellExperiment(list(counts=my.matrix))
#' gene.sets <- list( Tcell_signature = c("CD2","CD3E","CD3D"),
#'                  Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' my.sce <- ScoreSignatures_UCell_sce(my.sce, features=gene.sets)
#' head(t(assay(altExp(my.sce,"UCell"))))
#' ## End (Not run)
#' 
#' @export
ScoreSignatures_UCell_sce <- function(sce, assay="counts", features, ...) {
  
  if (class(sce) != "SingleCellExperiment") {
    stop("Provided object is not of class SingleCellExperiment.")
  }  
  if (!assay %in% names(assays(sce))) {
    stop(sprintf("Assay %s not found in sce object.", assay))
  }
  
  #Calculate UCell scores
  uscores <- ScoreSignatures_UCell(assay(sce, assay), features = features, ...)
  #Add scores to sce object
  altExp(sce, "UCell") <- SummarizedExperiment(assays = list("UCell" = t(uscores)))
  
  return(sce)
}
