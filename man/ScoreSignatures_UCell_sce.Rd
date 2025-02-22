% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScoreSignatures_UCell_sce.R
\name{ScoreSignatures_UCell_sce}
\alias{ScoreSignatures_UCell_sce}
\title{Calculate UCell scores on SingleCellExperiment objects}
\usage{
ScoreSignatures_UCell_sce(sce, assay = "counts", features, ...)
}
\arguments{
\item{sce}{A SingleCellExperiment object}

\item{assay}{The sce object assay where the data is to be found}

\item{features}{A list of signatures, for example: \code{list( Tcell_signature = c("CD2","CD3E","CD3D"), Myeloid_signature = c("SPI1","FCER1G","CSF1R"))}
You can also specify positive and negative gene sets by adding a + or - sign to genes in the signature; see an example below}

\item{...}{Additional parameters to be passed on to \code{\link{ScoreSignatures_UCell}}}
}
\value{
the SingleCellExperiment, with UCell scores added to altExp field
}
\description{
A wrapper for the function \code{\link{ScoreSignatures_UCell}} that interacts directly with SingleCellExperiment (sce) objects. (This could be an S4 generic)
}
\examples{
## Not run:
library(UCell)
my.matrix <- UCell::sample.matrix
my.sce <- SingleCellExperiment(list(counts=my.matrix))
gene.sets <- list( Tcell_signature = c("CD2","CD3E","CD3D"),
                 Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
my.sce <- ScoreSignatures_UCell_sce(my.sce, features=gene.sets)
head(t(assay(altExp(my.sce,"UCell"))))
## End (Not run)

}
