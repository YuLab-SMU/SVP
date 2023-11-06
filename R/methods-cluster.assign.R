#' @title clusting and assign the label for each feature(especifily the gene sets).
#' @rdname cluster.assign-method
#' @param data A \linkS4class{SVPExperiment}, which has run \code{sc.rwr} or \code{detect.svp}, or 
#' a \linkS4class{SingleCellExperiment} which was extracted from \linkS4class{SVPExperiment} using
#' \code{gsvaExp} function.
#' @param assay.type which expressed data to be pulled to run, default is \code{affi.score}.
#' @param ncluster integer the number cluster for each feature, default is 2.
#' @param random.seed integer random number to be reproducted.
#' @param BPPARAM A BiocParallelParam object specifying whether the built of KNN should be parallelized
#' default is \code{SerialParam()}, meaning no parallel. You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=T)}
#' to parallel it, the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param ... dot parameters
#' @export
setGeneric('cluster.assign', 
  function(
    data, 
    assay.type = 'affi.score',
    ncluster = 2,
    random.seed = 1024,
    BPPARAM = SerialParam(),
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
        ncluster = 2,
	random.seed = 1024,
        BPPARAM = SerialParam(),
	gsvaexp = NULL,
	gsvaexp.assay.type = NULL,
        ...
    ){
  if (is.null(assay.type)){
      assay.type <- 1
  }

  x <- assay(data, assay.type)
    
  y <- .internal.assign.cluster(x, ncluster, random.seed, BPPARAM, ...)

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
        ncluster = 2,
        random.seed = 1024,
        BPPARAM = SerialParam(),
	gsvaexp = NULL,
	gsvaexp.assay.type = NULL,
        ...
    ){
    
    if (!is.null(gsvaexp)){
       cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to clusting and assign.")
       da2 <- gsvaExp(data, gsvaexp)
       da2 <- cluster.assign(da2, assay.type = gsvaexp.assay.type, ncluster, random.seed, BPPARAM, ...)
       gsvaExp(data, gsvaexp) <- da2
    }else{
       data <- callNextMethod()
    }
    return(data)
    
})

#' @importFrom BiocParallel bplapply
#' @importFrom withr with_seed
#' @importFrom stats kmeans
.internal.assign.cluster <- function(da, k = 2, random.seed = 1024, BPPARAM = SerialParam(), ...){
    res <- with_seed(random.seed, 
               bplapply(seq(nrow(da)), 
                   function(i){
                     xx <- stats::kmeans(da[i, ], k, ...)

                     cluster <- xx$cluster
                     ind <- rank(xx$centers)

                     ind[cluster] - 1
                   },
                   BPPARAM = BPPARAM)
    )
    res <- do.call(rbind, res)
    rownames(res) <- rownames(da)
    colnames(res) <- colnames(da)
    res <- Matrix::Matrix(res, sparse = TRUE)
    return(res)
}

