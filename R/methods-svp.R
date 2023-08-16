#' @title svp
#' @rdname svp-methods
#' @param data object
#' @param gset.idx.list gene set list
#' @param .assay the assay name of \code{SingleCellExperiment}
#' @param gsva.method the method of gsva
#' @param gsva.kcdf the method of kcdf
#' @param gsva.abs.ranking logical whether rank the abs
#' @param gsva.min.sz integer the minimum gene set number
#' @param gsva.max.sz integer the maximum gene set number
#' @param gsva.parallel.sz the thread number
#' @param gsva.mx.diff logical 
#' @param gsva.tau exponent defining the weight of the tail in the random walk
#' @param gsva.ssgsea.norm logical whether normalizing when \code{gsva.method = 'ssgsea'}
#' @param gsva.verbose logical whether print the message of process in gsva
#' @param gsva.BPPARAM An object of class \code{BiocParallelParam} specifiying parameters
#' @param ... additional parameters
#' @importFrom BiocParallel SerialParam
#' @export
setGeneric('svp', function(data, 
                           gset.idx.list, 
                           .assay, 
                           gsva.method = c('gsva', 'ssgsea', 'zscore', 'plage'), 
                           gsva.kcdf = c('Gaussian', 'Poisson', 'none'), 
                           gsva.abs.ranking = FALSE,
                           gsva.min.sz = 1,
                           gsva.max.sz = Inf,
                           gsva.parallel.sz = 1L,
                           gsva.mx.diff = TRUE,
                           gsva.tau = switch(gsva.method, gsva=1, ssgsea = 0.25, NA),
                           gsva.ssgsea.norm = TRUE,
                           gsva.verbose = TRUE,
                           gsva.BPPARAM = SerialParam(progressbar = gsva.verbose),
                           ...)
    standardGeneric('svp')
)

#' @rdname svp-methods
#' @aliases svp,SingleCellExperiment,GeneSetCollection
#' @exportMethod svp
setMethod('svp', signature(data = 'SingleCellExperiment', gset.idx.list = 'GeneSetCollection'),
          function(
           data, 
           gset.idx.list,
           .assay,
           gsva.method = c('gsva', 'ssgsea', 'zscore', 'plage'),
           gsva.kcdf = c('Gaussian', 'Poisson', 'none'),
           gsva.abs.ranking = FALSE,
           gsva.min.sz = 1,
           gsva.max.sz = Inf,
           gsva.parallel.sz = 1L,
           gsva.mx.diff = TRUE,
           gsva.tau = switch(gsva.method, gsva=1, ssgsea = 0.25, NA),
           gsva.ssgsea.norm = TRUE,
           gsva.verbose = TRUE,
           gsva.BPPARAM = SerialParam(progressbar = gsva.verbose),
           ...){
           gsva.method <- match.arg(gsva.method)
           gsva.kcdf <- match.arg(gsva.kcdf)
           if (missing(.assay)).assay <- 1L
           x <- assay(data, .assay)
           x <- gsva(x, gset.idx.list, 
                       method = gsva.method, 
                       kcdf = gsva.kcdf,
                       abs.ranking = gsva.abs.ranking,
                       min.sz = gsva.min.sz,
                       max.sz = gsva.max.sz,
                       parallel.sz = gsva.parallel.sz,
                       mx.diff = gsva.mx.diff,
                       tau = gsva.tau,
                       ssgsea.norm = gsva.ssgsea.norm,
                       BPPARAM = gsva.BPPARAM
                  )
           x <- SingleCellExperiment(x)
           exp.nms <- paste0(ifelse(is.character(.assay), .assay, assayNames(data)[.assay]), ".", gsva.method)
           assayNames(x) <- exp.nms
           data <- .sce_to_svpe(data)
           gsvaExp(data, exp.nms) <- x
           return(data)
})

#' @rdname svp-methods
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @aliases svp,SingleCellExperiment,list
#' @exportMethod svp
setMethod('svp', signature(data = 'SingleCellExperiment', gset.idx.list = 'list'),
          function(
           data,
           gset.idx.list,
           .assay,
           gsva.method = c('gsva', 'ssgsea', 'zscore', 'plage'),
           gsva.kcdf = c('Gaussian', 'Poisson', 'none'),
           gsva.abs.ranking = FALSE,
           gsva.min.sz = 1,
           gsva.max.sz = Inf,
           gsva.parallel.sz = 1L,
           gsva.mx.diff = TRUE,
           gsva.tau = switch(gsva.method, gsva=1, ssgsea = 0.25, NA),
           gsva.ssgsea.norm = TRUE,
           gsva.verbose = TRUE,
           gsva.BPPARAM = SerialParam(progressbar = gsva.verbose),
           ...){
           gsva.method <- match.arg(gsva.method)
           gsva.kcdf <- match.arg(gsva.kcdf)
           if (missing(.assay)) .assay <- 1L
           x <- assay(data, .assay)
           x <- gsva(x, gset.idx.list,
                       annotation = .assay,
                       method = gsva.method,
                       kcdf = gsva.kcdf,
                       abs.ranking = gsva.abs.ranking,
                       min.sz = gsva.min.sz,
                       max.sz = gsva.max.sz,
                       parallel.sz = gsva.parallel.sz,
                       mx.diff = gsva.mx.diff,
                       tau = gsva.tau,
                       ssgsea.norm = gsva.ssgsea.norm,
                       BPPARAM = gsva.BPPARAM
                  ) |> suppressWarnings()
           x <- SingleCellExperiment(x)
           exp.nms <- paste0(ifelse(is.character(.assay), .assay, assayNames(data)[.assay]), ".", gsva.method)
           assayNames(x) <- exp.nms
           data <- .sce_to_svpe(data)
           gsvaExp(data, exp.nms) <- x
           return(data) 
})
