#' Detecting the spatially or single cell variable features with kullback–leibler divergence of 
#' 2D weighted kernel density estimation
#' @description
#' To resolve the sparsity of single cell or spatial omics data, we use kernel function smoothing cell 
#' density weighted by the gene expression in a low-dimensional space or physical space. This method had 
#' reported that it can better represent the gene expression, it can also recover the signal from cells that 
#' are more likely to express a gene based on their neighbouring cells (first reference). Next, we use
#' kullback-leibler divergence to detect the signal genes in a low-dimensional space (\code{UMAP} or \code{TSNE}
#' for single cell omics data) or a physical space (for spatial omics data). See details to learn more.
#'
#' @rdname runKldSVG-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param assay.type which expressed data to be pulled to run, default is \code{logcounts}.
#' @param reduction.used character used as spatial coordinates to detect SVG, default is \code{UMAP}, 
#' if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates.
#' @param sample_id character the sample(s) in the \linkS4class{SpatialExperiment} object whose cells/spots to use.
#' Can be \code{all} to compute metric for all samples; the metric is computed separately for each sample.
#' default is \code{"all"}.
#' @param grid.n numeric number of grid points in the two directions to estimate 2D weighted kernel 
#' density, default is 100.
#' @param permutation numeric the number of permutation for each single feature to detect the 
#' signicantly spatially or single cell variable features, default is 100.
#' @param p.adjust.method character the method to adjust the pvalue of the result, default is \code{BY}.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param action character control the type of output, if \code{action='add'}, the result of identification
#' will add the original object, if \code{action = 'get'}, the result will return a \linkS4class{SimpleList},
#' if \code{action = 'only'}, the result will return a \linkS4class{DataFrame} by merging the result of all
#' sample, default is \code{add}.
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
#'   \item \code{padj} the adjusted pvalue based on the speficied \code{p.adjust.method}, default is \code{BY}.
#'   \item \code{rank} the order of significant spatial variable features based on \code{padj} and \code{sp.kld}.
#' }
#' 
#' The kernel density estimation for each features in each cells is done in the following way (first reference article):
#'
#'   \eqn{f_{h}(x) = 1/n \sum_{i=1}^n W_{i} * K_{h}(x - X_{i})}
#'
#' Where \eqn{W_{i}} is the value of feature (such as gene expression or gene set score). \eqn{X_{i}} is the embeddings (two 
#' dimenstion coordinates of \code{UMAP} or \code{TSNE} or the physical space for spatial omics data) of the cell \eqn{i}.
#' \eqn{h} is a smoothing parameter corresponding to the bandwidth matrix, default is the implemention of \code{ks} package.
#' \eqn{K(x)} is a gaussian kernel function. \eqn{x} is the a reference point in the embedding space defined by the grid size 
#' used for the computation to weight the distances of nearby cells. \eqn{K_{h}(x—X_{i})}works as a weight for \eqn{W_{i}} to 
#' smooth the feature value based on neighbouring cells at a \code{UMAP} or \code{TSNE} or physical space.
#'
#' The Kullback-Leibler divergence for each features is calculated in the following way:
#' 
#'   \eqn{D_{KL}(G) = \sum_{x \in X} P(x) * \log(P(x) / Q(x))}
#' 
#'  Where \eqn{P(x)} is the kernel density value of a feature at the space \eqn{X}. and \eqn{Q(x)} is the kernel density value of no spatially
#'  variabilty reference feature at the space \eqn{X}. The smaller kullback-leibler divergence (\eqn{D_{KL}(G)}) show that the distribution of 
#'  features is more like the no spatially variabilty reference feature at th space \eqn{X}. So we randomly shuffle the position of each feature 
#'  and calculate Kullback-Leibler divergence, next we use the normal distribution to calculate the pvalue with the actual Kullback-Leibler 
#'  divergence, and the average value and standard deviation value of random Kullback-Leibler divergence, since the random Kullback-Leibler 
#'  divergence for each feature is normally distributed in the following:
#'  
#'   \eqn{X \sim \mathcal{N}(\mu,\,\sigma^{2})}
#' 
#'  where \eqn{\mu} is the average value of random Kullback-Leibler divergence, and \eqn{\sigma} is standard deviation.
#'
#' @references
#'
#' 1. Jose Alquicira-Hernandez, Joseph E Powell, Nebulosa recovers single-cell gene expression signals by kernel density estimation.
#'    Bioinformatics, 37, 2485–2487(2021), https://doi.org/10.1093/bioinformatics/btab003.
#' 
#' 2. Vandenbon, A., Diez, D. A clustering-independent method for finding differentially expressed genes in single-cell transcriptome 
#'    data. Nat Commun, 11, 4318 (2020). https://doi.org/10.1038/s41467-020-17900-3
#' 
#' 3. https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
#'
#' @seealso [`runSGSA`] to calculate the activity score of gene sets, [`runLISA`] to explore the hotspot for 
#' specified features in the spatial space.
#' @export
#' @author Shuangbin Xu
#' @examples
#' # This example dataset is extracted from the
#' # result of runSGSA with gsvaExp(svpe).
#' data(hpda_spe_cell_dec)
#' 
#' hpda_spe_cell_dec <-
#'     hpda_spe_cell_dec |>
#'     runKldSVG(
#'       assay.type = 'affi.score'
#'     )
#' 
#' # The result can be extracted svDf()
#' hpda_spe_cell_dec |> svDf() |> data.frame() |> dplyr::arrange(rank)
#' # the Acinar cells, Cancer clone A, Cancer clone B etc have
#' # significant spatial variable.
#' # Then we can use pred.feature.mode to predict the activity
#' # mode in spatial domain.
setGeneric('runKldSVG',
  function(
    data,
    assay.type = 'logcounts',
    reduction.used = c('UMAP', 'TSNE'),
    sample_id = 'all',
    grid.n = 100,
    permutation = 100,
    p.adjust.method = "BY",
    verbose = TRUE,
    action = c('add', 'only', 'get'),
    random.seed = 1024,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  )
  standardGeneric('runKldSVG')
)

#' @importFrom S4Vectors SimpleList
#' @rdname runKldSVG-method
#' @aliases runKldSVG,SingleCellExperiment
#' @export runKldSVG
setMethod('runKldSVG', 'SingleCellExperiment',
  function(
    data,
    assay.type = 'logcounts',
    reduction.used = c('UMAP', 'TSNE'),
    sample_id = 'all',
    grid.n = 100,
    permutation = 100,
    p.adjust.method = "BY",
    verbose = TRUE,
    action = c('add', 'only', 'get'),
    random.seed = 1024,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  ){
  action <- match.arg(action)
  if (is.null(assay.type)){
    assay.type <- assayNames(data)[1]
  }

  sample_id <- .check_sample_id(data, sample_id)

  x <- assay(data, assay.type)
  
  coords <- .check_coords(data, reduction.used, prefix="")

  tic()
  if (verbose){
      cli::cli_inform("Identifying the spatially variable gene sets (pathway) or genes based on
                       Kullback-Leibler divergence of 2D Weighted Kernel Density ...")

  }
  res.sv <- lapply(sample_id, function(sid){
              if (sid == ".ALLCELL"){
                 ind <- seq(ncol(x))
              }else{
                 ind <- colData(data)$sample_id == sid
              }
              xi <- x[, ind, drop=FALSE]
              xi <- xi[DelayedMatrixStats::rowVars(xi)!=0,]
              coordsi <- coords[ind,, drop=FALSE]
              grid.ni <- if(length(grid.n) > 1){grid.n[names(grid.n)==sid]}else{grid.n[1]}
              res.sv <- .identify.svg(
                                xi,
                                coords = coordsi,
                                n = grid.ni,
                                permutation = permutation,
                                p.adjust.method = p.adjust.method,
                                random.seed = random.seed,
                                ...)
            })
  toc()

  names(res.sv) <- sample_id

  if (action == 'get'){
      res.sv <- SimpleList(res.sv)
      return(res.sv)
  }

  res.sv <- .tidy_sv_result(res.sv)
  if (action == 'only'){
      return(res.sv)
  }
    
  if (verbose){
      cli::cli_inform(c("The result is added to the input object, which can be extracted using",
                       "`svDf()` with type='kld' for `SingleCellExperiment` or `SpatialExperiment`.", 
                       "If input object is `SVPExperiment`, and `gsvaexp` is specified, the result",
                       "can be extracted by `gsvaExp()` (return a `SingleCellExperiment` or ",
                       "`SpatialExperiment`, then also use `svDf()` to extract."))
  }

  data <- .add.int.rowdata(sce = data,
                        getfun = svDfs,
                        setfun1 = `svDfs<-`,
                        setfun2 = `svDf<-`,
                        namestr = 'sv.kld',
                        val = res.sv)

  return(data)
})

#' @rdname runKldSVG-method          
#' @aliases runKldSVG,SVPExperiment          
#' @export runKldSVG
setMethod('runKldSVG', 'SVPExperiment',
  function(
    data,
    assay.type = 'logcounts',
    reduction.used = c('UMAP', 'TSNE'),
    sample_id = 'all',
    grid.n = 100,
    permutation = 100,
    p.adjust.method = "BY",
    verbose = TRUE,
    action = c('add', 'only', 'get'),
    random.seed = 1024,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...){
    
    if (!is.null(gsvaexp)){
       if (verbose){
          cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to detect 'svg'.")
       }        
       da2 <- gsvaExp(data, gsvaexp, withSpatialCoords = TRUE, withReducedDim = TRUE, withColData = FALSE, withImgData=FALSE)
       da2 <- runKldSVG(da2, 
                     gsvaexp.assay.type, 
                     reduction.used,
                     sample_id, 
                     grid.n, 
                     permutation, 
                     p.adjust.method, 
                     verbose,
                     action, 
                     random.seed, 
                     ...)
       gsvaExp(data, gsvaexp) <- da2
    }else{
       data <- callNextMethod()
    }
    return(data)  
  }
)


#' Detecting the spatially or single cell variable features with Moran's I or Geary's C 
#' @description
#' This function use Moran's I, Geary's C or global G test to detect the signal genes in 
#' a low-dimensional space (\code{UMAP} or \code{TSNE} for single cell omics data) or 
#' a physical space (for spatial omics data).
#'
#' @rdname runDetectSVG-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param assay.type which expressed data to be pulled to run, default is \code{logcounts}.
#' @param method character one of \code{'moransi'}, \code{"gearysc"} or \code{"getisord"}, default is \code{'moransi'}.
#' @param weight object, which can be \code{nb}, \code{listw} or \code{Graph} object, default is NULL,
#' meaning the spatail neighbours weights will be calculated using the \code{weight.method}.
#' if the \code{data} contains multiple samples, and the \code{sample_id} is specified, it should be 
#' provided as a list object with names (using \code{sample_id}).
#' @param weight.method character the method to build the spatial neighbours weights, default 
#' is \code{knn} (k nearest neighbours). Other method, which requires coord matrix as input and returns
#' \code{nb}, \code{listw} or \code{Graph} object, also is avaiable, such as \code{'tri2nb'}, \code{"knearneigh"},
#' \code{'dnearneigh'}, \code{"gabrielneigh"}, \code{"relativeneigh"}, which are from \code{spdep} package.
#' default is \code{knn}, if it is \code{"none"}, meaning the distance weight of each spot is used to
#' the weight.
#' @param sample_id character the sample(s) in the \linkS4class{SpatialExperiment} object whose cells/spots to use. 
#' Can be \code{all} to compute metric for all samples; the metric is computed separately for each sample. 
#' default is \code{"all"}.
#' @param reduction.used character used as spatial coordinates to detect SVG, default is \code{UMAP},
#' if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates.
#' @param permutation integer the number to permutation test for the calculation of Moran's I, default
#' is NULL. Because we do not recommend using this parameter, as the permutation test is too slow.
#' @param p.adjust.method character the method to adjust the pvalue of the result, default is \code{BY}.
#' @param random.seed numeric random seed number to repeatability, default is 1024.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param action character control the type of output, if \code{action='add'}, the result of identification 
#' will add the original object, if \code{action = 'get'}, the result will return a \linkS4class{SimpleList},
#' if \code{action = 'only'}, the result will return a \linkS4class{DataFrame} by merging the result of all 
#' sample, default is \code{add}.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param ... additional parameters
#' @return a \linkS4class{SVPExperiment} or a \linkS4class{SingleCellExperiment}, see \code{action} parameter details.
#' @references
#' 1. Bivand, R.S., Wong, D.W.S. Comparing implementations of global and local indicators of spatial association. TEST 27, 
#'    716–748 (2018). https://doi.org/10.1007/s11749-018-0599-x
#' @export
#' @author Shuangbin Xu
#' @seealso [`runLISA`] to explore the hotspot for specified features in the spatial space.
#' @examples
#' # This example dataset is extracted from the
#' # result of runSGSA with gsvaExp(svpe).
#' data(hpda_spe_cell_dec)
#'
#' # using Moran's I test
#' ######################
#' hpda_spe_cell_dec <-
#'     hpda_spe_cell_dec |>
#'     runDetectSVG(
#'       assay.type = 'affi.score',
#'       method = 'moransi'
#'     )
#' # The result also is saved in the svDfs in the SVPExample object
#' # which can be extrated with svDf 
#' svDfs(hpda_spe_cell_dec)
#'
#' hpda_spe_cell_dec |> svDf("sv.moransi") |> data.frame() |> dplyr::arrange(rank)
#'
#' # using Geary's C test 
#' #######################
#' hpda_spe_cell_dec <-
#'     hpda_spe_cell_dec |>
#'     runDetectSVG(assay.type ='affi.score', method = 'gearysc')
#' 
#' svDfs(hpda_spe_cell_dec)
#' 
#' hpda_spe_cell_dec |> svDf("sv.gearysc") |> data.frame() |> dplyr::arrange(rank)
#'
#' # using Global G test (Getis-Ord)
#' #################################
#' hpda_spe_cell_dec <- hpda_spe_cell_dec |>
#'     runDetectSVG(assay.type = 1, method = 'getisord')
#'
#' svDfs(hpda_spe_cell_dec)
#' 
#' hpda_spe_cell_dec |> svDf(3) |> data.frame() |> dplyr::arrange(rank)
setGeneric("runDetectSVG", function(
    data,
    assay.type = 'logcounts',
    method = c("moransi", "gearysc", "getisord"),
    weight = NULL,
    weight.method = c("knn", "tri2nb", "none"),
    sample_id = "all",
    reduction.used = c('UMAP', 'TSNE'),
    permutation = NULL,
    p.adjust.method = "BH",
    random.seed = 1024,
    verbose = TRUE,
    action = c('add', 'only', 'get'),
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  )
  standardGeneric('runDetectSVG')
)

#' @rdname runDetectSVG-method
#' @aliases runDetectSVG,SingleCellExperiment
#' @export runDetectSVG
setMethod('runDetectSVG', 'SingleCellExperiment',
  function(
    data,
    assay.type = 'logcounts',
    method = c("moransi", "gearysc", "getisord"),
    weight = NULL,
    weight.method = c("knn", "tri2nb", "none"),
    sample_id = "all",
    reduction.used = c('UMAP', 'TSNE'),
    permutation = NULL,
    p.adjust.method = "BH",
    random.seed = 1024,
    verbose = TRUE,
    action = c('add', 'only', 'get'),
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...
  ){
  method <- match.arg(method)
  action <- match.arg(action)
  if (is.null(assay.type)){
    assay.type <- assayNames(data)[1]
  }else if (is.numeric(assay.type)){
    assay.type <- assayNames(data)[assay.type]
  }
  
  x <- assay(data, assay.type)

  sample_id <- .check_sample_id(data, sample_id)

  coords <- .check_coords(data, reduction.used, weight) 
  tic()
  if (verbose){
      cli::cli_inform(paste0("Identifying the spatially variable gene sets (pathway) based on ", method))
  }
  
  res.sv <- lapply(sample_id, function(sid){
      if (sid == ".ALLCELL"){
          ind <- seq(ncol(x))
      }else{
          ind <- colData(data)$sample_id == sid
      }
      xi <- x[, ind, drop=FALSE]
      xi <- xi[DelayedMatrixStats::rowVars(xi)!=0,]
      coordsi <- if(!is.null(coords)){coords[ind, , drop=FALSE]}else{NULL}
      weighti <- if(inherits(weight, "list")){weight[names(weight) == sid]}else{weight}
      res <- .identify.svg.by.autocorrelation(
                        xi,
                        coords = coordsi,
                        weight = weighti,
                        weight.method = weight.method,
                        method = method,
                        permutation = permutation,
                        p.adjust.method = p.adjust.method,
                        random.seed = random.seed,
                        ...)
      return(res)
  })

  toc()

  names(res.sv) <- sample_id

  if (action == 'get'){
      res.sv <- SimpleList(res.sv)
      return(res.sv)
  }

  res.sv <- .tidy_sv_result(res.sv)
  if (action == 'only'){
      return(res.sv)
  }
  nmstr <- paste0("sv.", method)
  if (verbose){
      cli::cli_inform(c("The result is added to the input object, which can be extracted using",
                       paste0("`svDf()` with type='", nmstr, "' for `SingleCellExperiment` or
                       `SpatialExperiment`."), "If input object is `SVPExperiment`, and `gsvaexp` is specified", 
                       "the result can be extracted by `gsvaExp()` (return a `SingleCellExperiment`",
                       " or `SpatialExperiment`),then also using `svDf()` to extract."))
  }

  data <- .add.int.rowdata(sce = data,
                        getfun = svDfs,
                        setfun1 = `svDfs<-`,
                        setfun2 = `svDf<-`,
                        namestr = paste0('sv.', method),
                        val = res.sv)

  return(data)  
})

#' @rdname runDetectSVG-method
#' @aliases runDetectSVG,SVPExperiment
#' @export runDetectSVG
setMethod('runDetectSVG', 'SVPExperiment',
  function(
    data,
    assay.type = 'logcounts',
    method = c("moransi", "gearysc"),
    weight = NULL,
    weight.method = c("knn", "tri2nb", "none"),
    sample_id = 'all',
    reduction.used = c('UMAP', 'TSNE'),
    permutation = NULL,
    p.adjust.method = "BH",
    random.seed = 1024,
    verbose = TRUE,
    action = c('add', 'only', 'get'),
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    ...){

    if (!is.null(gsvaexp)){
       if (verbose){
          cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to detect 'svg'.")
       }
       da2 <- gsvaExp(data, gsvaexp, withSpatialCoords = TRUE, withReducedDim = TRUE, withColData = FALSE, withImgData = FALSE)
       da2 <- runDetectSVG(da2,
                     gsvaexp.assay.type,
                     method,
                     weight,
                     weight.method,
                     sample_id,
                     reduction.used,
                     permutation,
                     p.adjust.method,
                     random.seed,
                     verbose,
                     action,
                     ...)
       gsvaExp(data, gsvaexp) <- da2
    }else{
       data <- callNextMethod()
    }
    return(data)
  }
)

