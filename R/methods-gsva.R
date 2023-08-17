#' @importFrom GSEABase geneIds
#' @export
setMethod('gsva',signature(expr='dgTMatrix', gset.idx.list='GeneSetCollection'), 
          function(expr, 
                   gset.idx.list, 
                   method = c("gsva", "ssgsea", "zscore", "plage"), 
                   kcdf = c("Gaussian", "Poisson", "none"),
                   abs.ranking = FALSE,
                   min.sz = 1,
                   max.sz = Inf,
                   parallel.sz = 1L,
                   mx.diff = TRUE,
                   tau = switch(method, gsva=1, ssgsea=0.25, NA),
                   ssgsea.norm = TRUE,
                   verbose = TRUE,
                   BPPARAM = SerialParam(progressbar=verbose),
                   ...){
    gset.idx.list <- GSEABase::geneIds(gset.idx.list)
    method <- match.arg(method)
    params <- list(...)
    params$expr <- expr
    params$gset.idx.list <- gset.idx.list
    params$method <- method
    params$kcdf <- kcdf
    params$abs.ranking <- abs.ranking
    params$min.sz <- min.sz
    params$max.sz <- max.sz
    params$parallel.sz <- parallel.sz
    params$mx.diff <- mx.diff
    params$tau <- tau
    params$ssgsea.norm <- ssgsea.norm
    params$verbose <- verbose
    params$BPPARAM <- BPPARAM
    do.call(gsva, params)
})

#' @importFrom methods as
#' @export
setMethod('gsva', signature(expr = 'dgTMatrix', gset.idx.list = 'list'),
          function(expr, gset.idx.list, 
                   method = c("gsva", "ssgsea", "zscore", "plage"),
                   kcdf = c("Gaussian", "Poisson", "none"),
                   abs.ranking = FALSE,
                   min.sz = 1,
                   max.sz = Inf,
                   parallel.sz = 1L,
                   mx.diff = TRUE,
                   tau = switch(method, gsva=1, ssgsea=0.25, NA),
                   ssgsea.norm = TRUE,
                   verbose = TRUE,
                   BPPARAM = SerialParam(progressbar=verbose),
                   ...){
    expr <- as(expr, 'CsparseMatrix')
    params <- list(...)
    method <- match.arg(method)
    params$expr <- expr
    params$gset.idx.list <- gset.idx.list
    params$method <- method
    params$kcdf <- kcdf
    params$abs.ranking <- abs.ranking
    params$min.sz <- min.sz
    params$max.sz <- max.sz
    params$parallel.sz <- parallel.sz
    params$mx.diff <- mx.diff
    params$tau <- tau
    params$ssgsea.norm <- ssgsea.norm
    params$verbose <- verbose
    params$BPPARAM <- BPPARAM
    do.call(gsva, params)
})
