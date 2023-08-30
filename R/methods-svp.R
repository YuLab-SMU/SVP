#' @title svp
#' @rdname svp-methods
#' @param data object
#' @param gset.idx.list gene set list
#' @param assayName the assay name of \code{SingleCellExperiment}
#' @param threads the thread number, default is 1.
#' @param gsva.method the method of gsva
#' @param gsva.kcdf the method of kcdf
#' @param gsva.abs.ranking logical whether rank the abs
#' @param gsva.min.sz integer the minimum gene set number
#' @param gsva.max.sz integer the maximum gene set number
#' @param gsva.mx.diff logical 
#' @param gsva.tau exponent defining the weight of the tail in the random walk
#' @param gsva.ssgsea.norm logical whether normalizing when \code{gsva.method = 'ssgsea'}
#' @param gsva.verbose logical whether print the message of process in gsva
#' @param gsva.BPPARAM An object of class \code{BiocParallelParam} specifiying parameters
#' @param sv.X matrix default is NULL.
#' @param sv.n_neighbors integer default is 10.
#' @param sv.order character default is \code{AMMD}.
#' @param sv.verbose logical
#' @param sv.cor.method character, the method to calculate the correlation of between 
#' the features and the color value for the RGB channels of histology image.
#' @param sv.cor.padj.method character, the method of p.adjust for the correlation test
#' @param sv.cor.image.beta integer, the square pixels to calculate the color value for
#' the RGB channels of histology image, default 2, meaning 2 x 2 pixels square.
#' @param ... additional parameters
#' @importFrom BiocParallel SerialParam
#' @export
setGeneric('svp', function(data, 
                           gset.idx.list, 
                           assayName,
                           threads = 1L, 
                           gsva.method = c('gsva', 'ssgsea', 'zscore', 'plage'), 
                           gsva.kcdf = c('Gaussian', 'Poisson', 'none'), 
                           gsva.abs.ranking = FALSE,
                           gsva.min.sz = 1,
                           gsva.max.sz = Inf,
                           gsva.mx.diff = TRUE,
                           gsva.tau = switch(gsva.method, gsva=1, ssgsea = 0.25, NA),
                           gsva.ssgsea.norm = TRUE,
                           gsva.verbose = TRUE,
                           gsva.BPPARAM = SerialParam(progressbar = gsva.verbose),
                           sv.X = NULL,
                           sv.n_neighbors = 10,
                           sv.order = 'AMMD',
                           sv.verbose = FALSE,
                           sv.cor.method = c("pearson", "spearman", "kendall"),
                           sv.cor.padj.method = 'fdr',
                           sv.cor.image.beta = 2,
                           ...)
    standardGeneric('svp')
)


#' @rdname svp-methods
#' @importFrom SummarizedExperiment rowData<- assayNames<-
#' @aliases svp,SingleCellExperiment,GeneSetCollection
#' @exportMethod svp
setMethod('svp', signature(data = 'SingleCellExperiment', gset.idx.list = 'GeneSetCollection'),
          function(
           data, 
           gset.idx.list,
           assayName,
           threads = 1L,
           gsva.method = c('gsva', 'ssgsea', 'zscore', 'plage'),
           gsva.kcdf = c('Gaussian', 'Poisson', 'none'),
           gsva.abs.ranking = FALSE,
           gsva.min.sz = 1,
           gsva.max.sz = Inf,
           gsva.mx.diff = TRUE,
           gsva.tau = switch(gsva.method, gsva=1, ssgsea = 0.25, NA),
           gsva.ssgsea.norm = TRUE,
           gsva.verbose = TRUE,
           gsva.BPPARAM = SerialParam(progressbar = gsva.verbose),
           sv.X = NULL,
           sv.n_neighbors = 10,
           sv.order = 'AMMD',
           sv.verbose = FALSE,           
           sv.cor.method = c('pearson', 'spearman', 'kendall'),
           sv.cor.padj.method = 'fdr',
           sv.cor.image.beta = 2,
           ...){
           gsva.method <- match.arg(gsva.method)
           gsva.kcdf <- match.arg(gsva.kcdf)
           if (missing(assayName))assayName <- 1L
           x <- assay(data, assayName)
           x <- gsva(x, gset.idx.list, 
                       method = gsva.method, 
                       kcdf = gsva.kcdf,
                       abs.ranking = gsva.abs.ranking,
                       min.sz = gsva.min.sz,
                       max.sz = gsva.max.sz,
                       parallel.sz = threads,
                       mx.diff = gsva.mx.diff,
                       tau = gsva.tau,
                       ssgsea.norm = gsva.ssgsea.norm,
                       BPPARAM = gsva.BPPARAM
                  ) |> suppressWarnings()

           flag <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)
           if(flag){
               specoords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
               res.sv <- .run_sv(x, spatial_coords = specoords, sv.X = sv.X, sv.order = sv.order,
                                 sv.n_neighbors = sv.n_neighbors, 
                                 sv.verbose = sv.verbose, 
                                 threads, ...)
               flag2 <- .check_element_obj(data, key = 'imgData', basefun = int_metadata, namefun = names)
               if (flag2){
                   sv.cor.method <- match.arg(sv.cor.method)
                   img <- .extract_element_object(data, key = 'imgData', basefun = int_metadata, namefun = names)
                   res.cor <- .cal_cor(x, img = img, coords = specoords,
                                       cor.method = sv.cor.method,
                                       beta = sv.cor.image.beta,
                                       threads = threads,
                                       p.adjust.method = sv.cor.padj.method,
                                       ...)
                   if (!is.null(res.cor)){
                       res.sv <- cbind(res.sv, res.cor)
                   }
               } 
           }

           x <- SingleCellExperiment(x)

           if(flag){
               res.sv <- .tidy_res.sv(x, res.sv)
               rowData(x) <- res.sv
           }

           exp.nms <- paste0(ifelse(is.character(assayName), assayName, assayNames(data)[assayName]), ".", gsva.method)
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
           assayName,
           threads = 1L,
           gsva.method = c('gsva', 'ssgsea', 'zscore', 'plage'),
           gsva.kcdf = c('Gaussian', 'Poisson', 'none'),
           gsva.abs.ranking = FALSE,
           gsva.min.sz = 1,
           gsva.max.sz = Inf,
           gsva.mx.diff = TRUE,
           gsva.tau = switch(gsva.method, gsva=1, ssgsea = 0.25, NA),
           gsva.ssgsea.norm = TRUE,
           gsva.verbose = TRUE,
           gsva.BPPARAM = SerialParam(progressbar = gsva.verbose),
           sv.X = NULL,
           sv.n_neighbors = 10,
           sv.order = 'AMMD',
           sv.verbose = FALSE,
           sv.cor.method = c('pearson', 'spearman', 'kendall'),
           sv.cor.padj.method = 'fdr',
           sv.cor.image.beta = 2,
           ...){
           gsva.method <- match.arg(gsva.method)
           gsva.kcdf <- match.arg(gsva.kcdf)
           if (missing(assayName)) assayName <- 1L
           x <- assay(data, assayName)
           x <- gsva(x, gset.idx.list,
                       annotation = assayName,
                       method = gsva.method,
                       kcdf = gsva.kcdf,
                       abs.ranking = gsva.abs.ranking,
                       min.sz = gsva.min.sz,
                       max.sz = gsva.max.sz,
                       parallel.sz = threads,
                       mx.diff = gsva.mx.diff,
                       tau = gsva.tau,
                       ssgsea.norm = gsva.ssgsea.norm,
                       BPPARAM = gsva.BPPARAM
                  ) |> suppressWarnings()

           flag <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)
           if(flag){
               specoords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
               res.sv <- .run_sv(x, spatial_coords = specoords, sv.X = sv.X, sv.order = sv.order, 
                                 sv.n_neighbors = sv.n_neighbors, sv.verbose = sv.verbose, threads = threads, 
                                 ...)
               flag2 <- .check_element_obj(data, key = 'imgData', basefun = int_metadata, namefun = names)
               if (flag2){
                   sv.cor.method <- match.arg(sv.cor.method)
                   img <- .extract_element_object(data, key = 'imgData', basefun = int_metadata, namefun = names) 
                   res.cor <- .cal_cor(x, img = img, coords = specoords, 
                                       cor.method = sv.cor.method,
                                       beta = sv.cor.image.beta, 
                                       threads = threads, 
                                       p.adjust.method = sv.cor.padj.method, 
                                       ...)
                   if (!is.null(res.cor)){
                       res.sv <- cbind(res.sv, res.cor)
                   }
               }
           }
           x <- SingleCellExperiment(x)
           if(flag){
               res.sv <- .tidy_res.sv(rowData(x), res.sv)
               rowData(x) <- res.sv
           }
           exp.nms <- paste0(ifelse(is.character(assayName), assayName, assayNames(data)[assayName]), ".", gsva.method)
           assayNames(x) <- exp.nms
           data <- .sce_to_svpe(data)
           gsvaExp(data, exp.nms) <- x
           return(data) 
})

#' @importFrom methods callNextMethod
#' @rdname svp-methods
#' @aliases svp,SpatialExperiment,GeneSetCollection
#' @exportMethod svp
setMethod('svp', signature(data = 'SpatialExperiment', gset.idx.list = 'GeneSetCollection'),
          function(
           data,
           gset.idx.list,
           assayName,
           threads = 1L,
           gsva.method = c('gsva', 'ssgsea', 'zscore', 'plage'),
           gsva.kcdf = c('Gaussian', 'Poisson', 'none'),
           gsva.abs.ranking = FALSE,
           gsva.min.sz = 1,
           gsva.max.sz = Inf,
           gsva.mx.diff = TRUE,
           gsva.tau = switch(gsva.method, gsva=1, ssgsea = 0.25, NA),
           gsva.ssgsea.norm = TRUE,
           gsva.verbose = TRUE,
           gsva.BPPARAM = SerialParam(progressbar = gsva.verbose),
           sv.X = NULL,
           sv.n_neighbors = 10,
           sv.order = 'AMMD',
           sv.verbose = FALSE,
           sv.cor.method = c('pearson', 'spearman', 'kendall'),
           sv.cor.padj.method = 'fdr',
           sv.cor.image.beta = 2, 
           ...){
           callNextMethod()
})

#' @rdname svp-methods
#' @aliases svp,SpatialExperiment,list
#' @exportMethod svp
setMethod('svp', signature(data = 'SpatialExperiment', gset.idx.list = 'list'),
          function(
           data,
           gset.idx.list,
           assayName,
           threads = 1L,
           gsva.method = c('gsva', 'ssgsea', 'zscore', 'plage'),
           gsva.kcdf = c('Gaussian', 'Poisson', 'none'),
           gsva.abs.ranking = FALSE,
           gsva.min.sz = 1,
           gsva.max.sz = Inf,
           gsva.mx.diff = TRUE,
           gsva.tau = switch(gsva.method, gsva=1, ssgsea = 0.25, NA),
           gsva.ssgsea.norm = TRUE,
           gsva.verbose = TRUE,
           gsva.BPPARAM = SerialParam(progressbar = gsva.verbose),
           sv.X = NULL,
           sv.n_neighbors = 10,
           sv.order = 'AMMD',
           sv.verbose = FALSE,
           sv.cor.method = c('pearson', 'spearman', 'kendall'),
           sv.cor.padj.method = 'fdr',
           sv.cor.image.beta = 2,
           ...){
           callNextMethod()
})



