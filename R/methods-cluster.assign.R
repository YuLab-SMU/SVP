#' @title clusting and assign the label for each feature(especifily the gene sets).
#' @rdname cluster.assign-method
#' @param data A \linkS4class{SVPExperiment}, which has run \code{runSGSA} or \code{detect.svp}, or 
#' a \linkS4class{SingleCellExperiment} which was extracted from \linkS4class{SVPExperiment} using
#' \code{gsvaExp} function.
#' @param assay.type which expressed data to be pulled to run, default is \code{affi.score}.
#' @param assign whether assign the max affinity of gene set or pathway to the each cell, default is FALSE.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param ... dot parameters
#' @return if input is a \linkS4class{SVPExperiment}, output will be also a \linkS4class{SVPExperiment}, and the result assay 
#' was stored in assay of the specified \code{gsvaexp}, which is a \linkS4class{SingleCellExperiment}. If input is a
#' \linkS4class{SingleCellExperiment} (which is extracted from \linkS4class{SVPExperiment} using \code{gsvaExp()} funtion), 
#' output will be a \linkS4class{SingleCellExperiment}, the result can be extracted using \code{assay()} function.
#' @details 
#' when use \code{runSGSA} to calculated the gene set activity of cell, if \code{assign = TRUE} we will assign the max affinity of
#' gene set or pathway to the each cell. If \code{assign = FALSE}, the max affinity of gene set or pathway will be kept.
#' @seealso to calculate the activity score of gene sets or pathway: [`runSGSA`].
#' @export
#' @examples
#' library(SpatialExperiment)
#' # This example data was extracted from the
#' # result of runSGSA with gsvaExp function.
#' data(hpda_spe_cell_dec)
#' assays(hpda_spe_cell_dec)
#' hpda_spe_cell_dec <- hpda_spe_cell_dec |>
#'    cluster.assign()
#' hpda_spe_cell_dec 
setGeneric('cluster.assign', 
  function(
    data, 
    assay.type = 'affi.score',
    assign = FALSE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL, 
    ...
  )
  standardGeneric('cluster.assign')
)

#' @rdname cluster.assign-method
#' @aliases cluster.assign,SingleCellExperiment
#' @export cluster.assign
setMethod(
    'cluster.assign', 
    'SingleCellExperiment', 
    function(
        data, 
        assay.type = 'affi.score',
        assign = FALSE,
        gsvaexp = NULL,
        gsvaexp.assay.type = NULL,
        ...
    ){
  if (is.null(assay.type)){
      assay.type <- 1
  }

  x <- assay(data, assay.type)
    
  y <- .internal.assign.cluster(x, assign, ...)

  assay(data, "cluster.assign")  <- y
  return(data)
})

#' @rdname cluster.assign-method
#' @aliases cluster.assign,SVPExperiment
#' @export cluster.assign
setMethod(
    'cluster.assign',
    'SVPExperiment',
    function(
        data,
        assay.type = 'affi.score',
        assign = FALSE,
        gsvaexp = NULL,
        gsvaexp.assay.type = NULL,
        ...
    ){
    
    if (!is.null(gsvaexp)){
       cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to clusting and assign.")
       da2 <- gsvaExp(data, gsvaexp, withColData = FALSE, withSpatialCoords=FALSE, withImageData = FALSE)
       da2 <- cluster.assign(da2, assay.type = gsvaexp.assay.type, assign, ...)
       gsvaExp(data, gsvaexp) <- da2
    }else{
       data <- callNextMethod()
    }
    return(data)
    
})

#' @importFrom withr with_seed
.internal.assign.cluster <- function(da, assign = FALSE, ...){
    i <- apply(da, 2, function(x) which.max(x))
    res <- sparseMatrix(i=i, j=seq(ncol(da)), x = 1, dims=c(nrow(da), ncol(da)))
    if (!assign){
        res <- da * res
    }
    rownames(res) <- rownames(da)
    colnames(res) <- colnames(da)
    return(res)
}

