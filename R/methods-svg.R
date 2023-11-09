#' @title Detecting the spatially or single cell variable features with Kullback–Leibler
#' divergence of 2D Weighted Kernel Density
#' @rdname kldSVG-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param assay.type which expressed data to be pulled to run, default is \code{logcounts}.
#' @param sv.used.reduction character used as spatial coordinates to detect SVG, default is \code{UMAP}, 
#' if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates.
#' @param sv.grid.n numeric number of grid points in the two directions to estimate 2D weighted kernel 
#' density, default is 100.
#' @param sv.permutation numeric the number of permutation for each single feature to detect the 
#' signicantly spatially or single cell variable features, default is 100.
#' @param sv.p.adjust.method character the method to adjust the pvalue of the result, default is \code{bonferroni}.
#' @param sv.BPPARAM A BiocParallelParam object specifying whether the identification of SV features should be parallelized
#' default is \code{SerialParam()}, meaning no parallel. You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=T)}
#' to parallel it, the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}. default is \code{SerialParam()}.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param random.seed numeric random seed number to repeatability, default is 1024.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param ... additional parameters
#' @return a \linkS4class{SVPExperiment} or a \linkS4class{SingleCellExperiment}, see details.
#' @details
#' if input is a \linkS4class{SVPExperiment}, output will be also a \linkS4class{SVPExperiment}, the spatially variable gene sets 
#' result is stored in \code{svDfs} of the specified \code{gsvaexp}, which is a \linkS4class{SingleCellExperiment}. If input is 
#' a \linkS4class{SingleCellExperiment} (which is extracted from \linkS4class{SVPExperiment} using \code{gsvaExp()} funtion), output
#' will be also a \linkS4class{SingleCellExperiment}, the spatial variable gene sets result can be extracted using \code{svDf} function.
#' The result of \code{svDf} will return a matrix which has \code{sp.kld}, \code{boot.sp.kld.mean}, \code{boot.sp.kld.sd}, \code{pvalue},
#' \code{padj} and \code{rank}.
#' \itemize{
#'   \item \code{sp.kld} which is logarithms of Kullback–Leibler divergence, larger value meaning the greater the difference from the
#'      background distribution without spatial variability.
#'   \item \code{boot.sp.kld.mean} which is mean of logarithms of Kullback–Leibler divergence based on the permutation of each features.
#'   \item \code{boot.sp.kld.sd} which is standard deviation of logarithms of Kullback–Leibler divergence based on the permutation of
#'      each features.
#'   \item \code{pvalue} the pvalue is calculated using the real \code{sp.kld} and the permutation \code{boot.sp.kld.mean} and
#'      \code{boot.sp.kld.sd} based on the normal distribution.
#'   \item \code{padj} the adjusted pvalue based on the speficied \code{sv.p.adjust.method}, default is \code{BY}.
#'   \item \code{rank} the order of significant spatial variable features based on \code{padj} and \code{sp.kld}.
#' }
#' @seealso [`sc.rwr`] to calculate the activity score of gene sets.
#' @export
setGeneric('kldSVG',
  function(
    data,
    assay.type = 'logcounts',
    sv.used.reduction = c('UMAP', 'TSNE'),
    sv.grid.n = 100,
    sv.permutation = 100,
    sv.p.adjust.method = "bonferroni",
    sv.BPPARAM = SerialParam(),
    verbose = TRUE,
    random.seed = 1024,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  )
  standardGeneric('kldSVG')
)

#' @rdname kldSVG-method
#' @aliases kldSVG,SingleCellExperiment
#' @export kldSVG
setMethod('kldSVG', 'SingleCellExperiment',
  function(
    data,
    assay.type = 'logcounts',
    sv.used.reduction = c('UMAP', 'TSNE'),
    sv.grid.n = 100,
    sv.permutation = 100,
    sv.p.adjust.method = "bonferroni",
    sv.BPPARAM = SerialParam(),
    verbose = TRUE,
    random.seed = 1024,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  ){
  if (is.null(assay.type)){
    assay.type <- assayNames(data)[1]
  }

  x <- assay(data, assay.type)

  flag1 <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)

  flag2 <- any(sv.used.reduction %in% reducedDimNames(data))

  if((flag1 || flag2)){
      if (flag2){
          coords <- reducedDim(data, sv.used.reduction)
          coords <- coords[,c(1, 2)]
      }
      if (flag1){
          coords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
      }
      tic()
      if (verbose){
          cli::cli_inform("Identifying the spatially variable gene sets (pathway) based on
                           Kullback-Leibler divergence of 2D Weighted Kernel Density ...")

      }

      res.sv <- .identify.svg(
                        x,
                        coords = coords,
                        n = sv.grid.n,
                        permutation = sv.permutation,
                        p.adjust.method = sv.p.adjust.method,
                        BPPARAM = sv.BPPARAM,
                        random.seed = random.seed,
                        ...)

      data <- .add.int.rowdata(sce = data,
                            getfun = svDfs,
                            setfun1 = `svDfs<-`,
                            setfun2 = `svDf<-`,
                            namestr = 'sv.kld',
                            val = res.sv)
      toc()
  }else{
      cli::cli_abort("The {.cls {class(data)}} should have 'spatialCoords' or the reduction result of 'UMAP' or 'TSNE'.")
  }

  return(data)
})

#' @rdname kldSVG-method          
#' @aliases kldSVG,SVPExperiment          
#' @export kldSVG
setMethod('kldSVG', 'SVPExperiment',
  function(
    data,
    assay.type = 'logcounts',
    sv.used.reduction = c('UMAP', 'TSNE'),
    sv.grid.n = 100,
    sv.permutation = 100,
    sv.p.adjust.method = "bonferroni",
    sv.BPPARAM = SerialParam(),
    verbose = TRUE,
    random.seed = 1024,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...){
    
    if (!is.null(gsvaexp)){
       if (verbose){
          cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to detect 'svg'.")
       }        
       da2 <- gsvaExp(data, gsvaexp, withSpatialCoords = TRUE, withReducedDim = TRUE)
       da2 <- kldSVG(da2, 
                     gsvaexp.assay.type, 
                     sv.used.reduction, 
                     sv.grid.n, 
                     sv.permutation, 
                     sv.p.adjust.method, 
                     sv.BPPARAM, 
                     verbose, 
                     random.seed, 
                     ...)
       gsvaExp(data, gsvaexp, withSpatialCoords = FALSE, withReducedDim = FALSE) <- da2
    }else{
       data <- callNextMethod()
    }
    return(data)  
  }
)
