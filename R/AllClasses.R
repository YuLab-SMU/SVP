#' @title The SVPExperiment class
#' @docType class
#' @param ... passed to the \code{\link{SingleCellExperiment}} constructor to 
#' fill the slots of the base class.
#' @param gsvaExps list containing \linkS4class{SingleCellExperiment} object,
#' each of which should have the same number of columns as the output 
#' \linkS4class{SVPExperiment} object.
#' @importFrom methods setClass setClassUnion
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @return a \linkS4class{SVPExperiment} object
#' @aliases
#' coerce,SingleCellExperiment,SVPExperiment-method
#' @export
#' @author Shuangbin Xu
#' @examples
#' library(SingleCellExperiment) |> suppressPackageStartupMessages()
#' ncells <- 100
#' u <- matrix(rpois(20000, 5), ncol=ncells)
#' v <- log2(u + 1)
#' pca <- matrix(runif(ncells*5), ncells)
#' tsne <- matrix(rnorm(ncells*2), ncells)
#'
#' svpe <- SVPExperiment(assays=list(counts=u, logcounts=v),
#'     reducedDims=SimpleList(PCA=pca, tSNE=tsne))
#' svpe
#'
#' ## coercion from SingleCellExperiment
#' sce <- SingleCellExperiment(assays=list(counts=u, logcounts=v),
#' reducedDims=SimpleList(PCA=pca, tSNE=tsne))
#' svpe <- as(sce, 'SVPExperiment')
#' svpe
SVPExperiment <- function(..., gsvaExps = list()){
    svpe <- SingleCellExperiment(...)
    svpe <- .sce_to_svpe(svpe, gsvaExps = gsvaExps)
    return(svpe)
}


#' @export
#' @rdname SVPExperiment
setClass('SVPExperiment', contains = 'SingleCellExperiment')

setClass('SCEByColumn', slot = c(sce = 'SingleCellExperiment'))

setClassUnion("matrix_Or_NULL", c("matrix", "NULL"))

#' @exportMethod coerce
#' @importFrom methods as
setAs(from = "SingleCellExperiment", to = 'SVPExperiment', function(from){
    svpe <- .sce_to_svpe(from)
    return(svpe)
})
