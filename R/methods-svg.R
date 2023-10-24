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
#' @param sv.p.adjust.method character the method to adjust the pvalue of the result, default is \code{fdr}.
#' @param sv.BPPARAM A BiocParallelParam object specifying whether the identification of SV features should be parallelized
#' default is \code{SerialParam()}, meaning no parallel. You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=T)}
#' to parallel it, the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}. default is \code{SerialParam()}.
#' @param cells Vector specifying the subset of cells to be used for the calculation of the activaty score or identification
#' of SV features. This can be a character vector of cell names, an integer vector of column indices or a logical vector,
#' default is NULL, meaning all cells to be used for the calculation of the activaty score or identification of SV features.
#' @param features Vector specifying the subset of features to be used for the calculation of the activaty score or identification
#' of SV features. This can be a character vector of features names, an integer vector of row indices or a logical vector,
#' default is NULL, meaning all features to be used for the calculation of the activaty score or identification of SV features.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param random.seed numeric random seed number to repeatability, default is 1024.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param ... additional parameters
#' @export
setGeneric('kldSVG',
  function(
    data,
    assay.type = 'logcounts',
    sv.used.reduction = c('UMAP', 'TSNE'),
    sv.grid.n = 100,
    sv.permutation = 100,
    sv.p.adjust.method = "fdr",
    sv.BPPARAM = SerialParam(),
    cells = NULL,
    features = NULL,
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
    sv.p.adjust.method = "fdr",
    sv.BPPARAM = SerialParam(),
    cells = NULL,
    features = NULL,
    verbose = TRUE,
    random.seed = 1024,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  ){
  if (!assay.type %in% assayNames(data)){
      cli::cli_abort("the {.var assay.type} = {assay.type} is not present in the assays of {.cls {class(data)}}.")
  }

  if (!is.null(cells)){
    data <- data[cells,]
  }

  if (!is.null(features)){
    data <- data[,features]
  }

  x <- assay(data, assay.type)

  flag <- .check_element_obj(data, key = 'gsvaExps', basefun=int_colData, namefun = names)
  if (!is.null(gsvaexp) && flag){
    if (verbose){
       cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to detect 'svg'.")
    }
    data <- gsvaExp(data, gsvaexp, withSpatialCoords = TRUE, withReducedDim = TRUE)
    if (is.null(gsvaexp.assay.type)){gsvaexp.assay.type <- 1}
    x <- assay(data, gsvaexp.assay.type)
  }

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
    sv.p.adjust.method = "fdr",
    sv.BPPARAM = SerialParam(),
    cells = NULL,
    features = NULL,
    verbose = TRUE,
    random.seed = 1024,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...){
    
    x <- callNextMethod()
    
    if (!is.null(gsvaexp)){
       gsvaExp(data, gsvaexp, withSpatialCoords = FALSE, withReducedDim = FALSE) <- x
       return(data)
    }
    return(x)  
  }
)
