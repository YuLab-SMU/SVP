#' @title predict the cell signature according the gene sets or pathway activity score.
#' @rdname pred.cell.signature-method
#' @param data A \linkS4class{SVPExperiment}, which has run \code{runSGSA} or \code{detect.svp}, or 
#' a \linkS4class{SingleCellExperiment} which was extracted from \linkS4class{SVPExperiment} using
#' \code{gsvaExp} function.
#' @param assay.type which expressed data to be pulled to run, default is \code{affi.score}.
#' @param threshold numeric when the gene set activity score of cell less than the \code{threshold}, the cell
#' signature will be consider as 'unassigned', default is NULL, meaning will be calculated internally.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param pred.col.name character the column name in \code{colData} of the result, default is \code{pred.cell.sign}.
#' @param ... dot parameters
#' @return if input is a \linkS4class{SVPExperiment}, output will be also a \linkS4class{SVPExperiment}, and the result 
#' was stored at the \code{pred.col.name} column of \code{colData} in the specified \code{gsvaexp}, which is a 
#' \linkS4class{SingleCellExperiment}. If input is a \linkS4class{SingleCellExperiment} (which is extracted from 
#' \linkS4class{SVPExperiment} using \code{gsvaExp()} funtion), output will be a \linkS4class{SingleCellExperiment}, 
#' the result can be extracted using \code{colData()} function with specified column in default is \code{pred.cell.sign}.
#' @seealso to calculate the activity score of gene sets or pathway: [`runSGSA`], 
#' to keep the max gene set or pathway activity score of cell: [`cluster.assign`].
#' @export
setGeneric('pred.cell.signature', 
  function(
    data, 
    assay.type = 'affi.score',
    threshold = NULL,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL, 
    pred.col.name = 'pred.cell.sign',
    ...
  )
  standardGeneric('pred.cell.signature')
)

#' @rdname pred.cell.signature-method
#' @aliases pred.cell.signature,SingleCellExperiment
#' @export pred.cell.signature
setMethod(
    'pred.cell.signature', 
    'SingleCellExperiment', 
    function(
        data, 
        assay.type = 'affi.score',
        threshold = NULL,
        gsvaexp = NULL,
        gsvaexp.assay.type = NULL,
	pred.col.name = 'pred.cell.sign',
        ...
    ){
  if (is.null(assay.type)){
      assay.type <- 1
  }

  x <- assay(data, assay.type)
    
  colData(data)[[pred.col.name]] <- .internal.predict.cell.sign(x, threshold, ...)

  return(data)
})

#' @rdname pred.cell.signature-method
#' @aliases pred.cell.signature,SVPExperiment
#' @export pred.cell.signature
setMethod(
    'pred.cell.signature',
    'SVPExperiment',
    function(
        data,
        assay.type = 'affi.score',
        threshold = NULL,
        gsvaexp = NULL,
        gsvaexp.assay.type = NULL,
        pred.col.name = 'pred.cell.sign',
        ...
    ){
    
    if (!is.null(gsvaexp)){
       cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to predict the cell signature.")
       da2 <- gsvaExp(data, gsvaexp, withColData=FALSE, withSpatialCoords=FALSE, withImageData = FALSE)
       da2 <- pred.cell.signature(da2, assay.type = gsvaexp.assay.type, threshold = threshold, pred.col.name = pred.col.name, ...)
       gsvaExp(data, gsvaexp) <- da2
    }else{
       data <- callNextMethod()
    }
    return(data)
})

#' @importFrom BiocParallel bplapply
#' @importFrom withr with_seed
.internal.predict.cell.sign <- function(da, threshold = NULL, ...){
    pred <- rownames(da)[apply(da, 2, function(x) which.max(x))]
    if (is.null(threshold)){
        threshold <- 0
    }
    pred <- ifelse(apply(da, 2, max) > threshold, pred, 'unassigned')
    return(pred)
}

