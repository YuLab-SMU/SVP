#' @title Local indicators of spatial association analysis 
#' @description
#' This function use the local indicators of spatial association (LISA) to identify the hotspot
#' in the spatial space.
#' @rdname runLISA-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param features the feature name or index of data object, which are required. If \code{gsvaexp} is
#' provided and \code{data} is \linkS4class{SingleCellExperiment}, it should be the features from
#' `rownames(gsvaExp(data, gsvaexp))`.
#' @param assay.type which expressed data to be pulled to run, default is \code{logcounts}.
#' @param sample_id character the sample(s) in the \linkS4class{SpatialExperiment} object whose cells/spots to use.
#' Can be \code{all} to compute metric for all samples; the metric is computed separately for each sample.
#' default is \code{"all"}.
#' @param method character the method for the local spatial statistic, one of \code{'localG'}, \code{"localmoran"}, 
#' default is \code{'localG'}.
#' @param weight object, which can be \code{nb}, \code{listw} or \code{Graph} object, default is NULL,
#' meaning the spatail neighbours weights will be calculated using the \code{weight.method}.
#' if the \code{data} contains multiple samples, and the \code{sample_id} is specified, it should be
#' provided as a list object with names (using \code{sample_id}).
#' @param weight.method character the method to build the spatial neighbours weights, default
#' is \code{voronoi} (Voronoi tessellation). Other method, which requires coord matrix as input and returns
#' \code{nb}, \code{listw} or \code{Graph} object, also is avaiable, such as \code{"knearneigh"},
#' \code{'dnearneigh'}, \code{"gabrielneigh"}, \code{"relativeneigh"}, which are from \code{spdep} package.
#' default is \code{knn}, if it is \code{"none"}, meaning the distance weight of each spot is used to
#' the weight.
#' @param reduction.used character used as spatial coordinates to calculate the neighbours weights,
#' default is \code{NULL}, the result of reduction can be specified, such as \code{UMAP}, \code{TSNE}, \code{PCA}.
#' If it is specified, the weight neighbours matrix will be calculated using the result of specified reduction.
#' @param group.by character a specified category column names (for example the cluster column name) of
#' \code{colData(data)}, if it was specified, the adjacency weighted matrix will be built based on the principle
#' that spots or cells in the same category are adjacent, default is NULL.
#' @param cells the cell name or index of data object, default is NULL.
#' @param action character, which control the type of return result, default is \code{get}, which will return
#' a \linkS4class{SimpleList}.
#' @param alternative a character string specifying the alternative hypothesis, default is \code{two.sided}.
#' @param BPPARAM A BiocParallelParam object specifying whether perform the analysis parallelly using
#' \code{BiocParallel} default is \code{SerialParam()}, meaning no parallel.
#' You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=T)} to parallel it,
#' the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}. default is \code{SerialParam()}.
#' @param verbose logical whether print the help information, default is TRUE.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param gsvaexp.features character which is from the `rownames(gsvaExp(data, gsvaexp))`. If \code{gsvaexp} is
#' specified and \code{data} is \linkS4class{SVPExperiment}, it should be provided. Default is NULL.
#' @param ... additional parameters the parameters which are from the weight.method function.
#' @return if \code{action = 'get'} (in default), the SimpleList object (like list object) will be return,
#' if \code{action = 'only'}, the data.frame will be return. if \code{action = 'add'}, the result of LISA is
#' stored in the \code{localResults} column of \code{int_colData} (internal column metadata), which can be extracted
#' using [`LISAResult`]
#' @references
#' 1. Bivand, R.S., Wong, D.W.S. Comparing implementations of global and local indicators of spatial association. TEST 27,
#'    716–748 (2018). https://doi.org/10.1007/s11749-018-0599-x
#' @seealso [`runDetectSVG`] and [`runKldSVG`] to identify the spatial variable features.
#' @author Shuangbin Xu
#' @export
#' @examples
#' library(SpatialExperiment)
#' # This example data was extracted from the
#' # result of runSGSA with gsvaExp() function.
#' data(hpda_spe_cell_dec)
#' # using global spatial autocorrelation test to identify the spatial 
#' # variable features.
#' svres <- runDetectSVG(hpda_spe_cell_dec, assay.type = 'affi.score', 
#'            method = 'moransi', action = 'only') 
#' svres |> dplyr::arrange(rank) |> head()
#' # In this example, we found the `Cancer clone A` and `Cancer clone B`
#' # have significant spatial autocorrelation. Next, we use the `runLISA()`
#' # to explore the spatial hotspots for the features.
#' lisa.res12 <- hpda_spe_cell_dec |>
#'    runLISA(
#'      features = c(1, 2, 3), 
#'      assay.type = 'affi.score',
#'      weight.method = "knn",
#'      k = 10,
#'      action = 'get',
#'    )
#' lisa.res12
#' lisa.res12[['Acinar cells']] |> head()
#' lisa.res12[["Cancer clone A"]] |> head()
#' colData(hpda_spe_cell_dec)$`cluster.test.Cancer.A` <- lisa.res12[["Cancer clone A"]] |>
#' dplyr::pull(cluster.test)
#' colData(hpda_spe_cell_dec)$`cluster.test.Acinar` <- lisa.res12[["Acinar cells"]] |>
#' dplyr::pull(cluster.test)
#' colData(hpda_spe_cell_dec)$`cluster.test.Cancer.B` <- lisa.res12[["Cancer clone B"]] |>
#' dplyr::pull(cluster.test)
#' # Then using ggsc to visualize the result
#' \dontrun{
#'   library(ggplot2)
#'   library(ggsc)
#'   p1 <- sc_spatial(hpda_spe_cell_dec,
#'              features = rownames(hpda_spe_cell_dec),
#'              mapping = aes(x=x, y=y, color=cluster.test.Cancer.A),
#'              plot.pie = T,
#'              pie.radius.scale = .8,
#'              bg_circle_radius = 1.1,
#'              color=NA,
#'              linewidth=2
#'   ) +
#'   scale_color_manual(values=c("black", "white"))
#'   p1
#'   f1 <- sc_spatial(hpda_spe_cell_dec, features="Cancer clone A",
#'              mapping=aes(x = x, y = y),
#'              pointsize=10
#'   ) +
#'   geom_scattermore2(
#'     mapping = aes(bg_color=cluster.test.Cancer.A, subset=cluster.test.Cancer.A=="High"),
#'     bg_line_width = .16,
#'     gap_line_width = .14,
#'     pointsize = 7
#'   ) +
#'   scale_bg_color_manual(values=c('black'))
#'   f1
#'   f2 <- sc_spatial(hpda_spe_cell_dec, features="Acinar cells",
#'              mapping=aes(x=x,y=y),
#'              pointsize=10
#'   ) +
#'   geom_scattermore2(
#'     mapping = aes(bg_color=cluster.test.Acinar, subset=cluster.test.Acinar=="High"),
#'     bg_line_width = .16,
#'     gap_line_width = .14,
#'     pointsize = 7
#'   ) +
#'   scale_bg_color_manual(values=c('black'))
#'   f2
#'   f3 <- sc_spatial(hpda_spe_cell_dec, features="Cancer clone B",
#'              mapping=aes(x=x,y=y),
#'              pointsize=10
#'   ) +
#'   geom_scattermore2(
#'     mapping = aes(bg_color=cluster.test.Cancer.B, subset=cluster.test.Cancer.B=="High"),
#'     bg_line_width = .18,
#'     gap_line_width = .14,
#'     pointsize = 8
#'   ) +
#'   scale_bg_color_manual(values=c('black'))
#' }
setGeneric('runLISA',
  function(
    data,
    features,
    assay.type = 'logcounts',
    sample_id = 'all',
    method = c("localG", "localmoran"),
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    reduction.used = NULL,
    group.by = NULL,
    cells = NULL,
    action = c("get", "add", "only"),
    alternative = 'two.sided',
    BPPARAM = SerialParam(),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    ...
  )
  standardGeneric('runLISA')
)

#' @rdname runLISA-method
#' @aliases runLISA,SingleCellExperiment
#' @export runLISA
setMethod("runLISA", "SingleCellExperiment", function(
    data, 
    features, 
    assay.type = "logcounts",
    sample_id = 'all',
    method = c("localG", "localmoran"),
    weight = NULL, 
    weight.method = c("voronoi", "knn", "none"), 
    reduction.used = NULL,
    group.by = NULL,
    cells = NULL,
    action = c("get", "add", "only"),
    alternative = 'two.sided',
    BPPARAM = SerialParam(),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    ...
  ){
  
  action <- match.arg(action)
  weight.method <- match.arg(weight.method)
  method <- match.arg(method)
  sample_id <- .check_sample_id(data, sample_id)
  
  if (is.null(assay.type)){
      assay.type <- 1
  }

  x <- assay(data, assay.type)

  if (!is.null(cells)){
      x <- x[, cells, drop=FALSE]
  }

  x <- x[features, ,drop=FALSE]
  
  weight <- .check_weight(data, sample_id, weight, group.by)
  coords <- .check_coords(data, reduction.used, weight, weight.method)

  res <- lapply(sample_id, function(sid){
                  if (sid == ".ALLCELL"){
                      ind <- seq(ncol(x))
                  }else{
                      ind <- colData(data)$sample_id == sid
                  }
                  coordsi <- if(!is.null(coords)){coords[ind, ,drop=FALSE]}else{NULL}
                  weighti <- if(inherits(weight, 'list')){weight[names(weight) == sid]}else{weight}
                  xi <- x[, ind, drop=FALSE]
                  wm <- .obtain.weight(coordsi, weight = weighti, weight.method = weight.method, ...)
                  n <- nrow(wm)
                  wi <- rowSums(wm)
                  wi2 <- rowSums(wm^2)
                  if (any(wi == 0)){
                      cli::cli_warn("no-neighbour observations found in the spatial neighborhoods graph.")
                  }
                  res <- BiocParallel::bplapply(seq(nrow(xi)), function(i){ 
                                                .internal.runLISA(xi[i,], wm, wi, wi2, n, method, alternative)
                             }, BPPARAM = BPPARAM           
                         )
                  names(res) <- rownames(x)
                  return(res)
         })
  res <- .tidy_lisa_res(res)

  if (action == 'only'){
      res <- do.call("cbind", res)
      return(res)
  }

  if (action == 'add'){
      if (verbose){
          cli::cli_inform(c("The result is added to the input object, which can be extracted using",
                           paste0("`LISAResult()` with `type='", method, ".SVP'`, and a specified `features`.")))
      }
      data <- .add_LISA_res(data, res, method)
      return(data)
  }
  res <- S4Vectors::SimpleList(res)
  return(res)
})



#' @rdname runLISA-method
#' @aliases runLISA,SVPExperiment
#' @export runLISA
setMethod("runLISA", "SVPExperiment", function(
    data,
    features,
    assay.type = "logcounts",
    sample_id = 'all',
    method = c("localG", "localmoran"),
    weight = NULL,
    weight.method = c("voronoi", "knn", "none"),
    reduction.used = NULL,
    group.by = NULL,
    cells = NULL,
    action = c("get", "add", "only"),
    alternative = 'two.sided',
    BPPARAM = SerialParam(),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    ...
  ){
    method <- match.arg(method)
    action <- match.arg(action)
    if (!is.null(gsvaexp)){
        if (verbose){
            cli::cli_inform(c("The {.var gsvaexp} is specified, the specified {.var gsvaExp} will be used to perform LISA."))
        }
        if (missing(features) && is.null(gsvaexp.features)){
            cli::cli_abort("The {.var gsvaexp} is specified, the {.var features} or {.var gsvaexp.features} should be provided.")
        }
        if (missing(features)){
            features <- gsvaexp.features
        }else{
            features <- unique(features, gsvaexp.features)
        }
        da2 <- gsvaExp(data, gsvaexp, withSpatialCoords = TRUE, withReducedDim = TRUE, withColData = TRUE, withImgData = FALSE)
        da2 <- runLISA(da2,
                       features,
                       gsvaexp.assay.type,
                       sample_id,
                       method,
                       weight,
                       weight.method,
                       reduction.used,
                       group.by,
                       cells,
                       action,
                       alternative,
                       BPPARAM,       
                       verbose,
                      ...)
        if (action!='add'){
            return(da2)
        }
        if (verbose){
            cli::cli_inform(c("The `action = 'add'`, If you want to extract the results, you need to",
                              "use `gsvaExp(data, gsvaexp)` to extract the internal object firstly."))
        }
        gsvaExp(data, gsvaexp) <- da2
    }else{
        data <- callNextMethod()
    }
    return(data)
})
