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


#' @title predict the expression or activity score mode of the 
#' gene sets or pathway activity score or other features.
#' @description
#' this function will return a certain number of modes for each features using the estimation of 
#' the location of modes and antimodes with the expression or activity score.
#' @rdname pred.feature.mode-method
#' @param data A \linkS4class{SVPExperiment}, which has run \code{runSGSA} or \code{detect.svp}, or
#' a \linkS4class{SingleCellExperiment} which was extracted from \linkS4class{SVPExperiment} using
#' \code{gsvaExp} function.
#' @param assay.type which expressed data to be pulled to run, default is \code{affi.score}.
#' @param gsvaexp which gene set variation experiment will be pulled to run, this only work when \code{data} is a
#' \linkS4class{SVPExperiment}, default is NULL.
#' @param gsvaexp.assay.type which assay data in the specified \code{gsvaexp} will be used to run, default is NULL.
#' @param mod0 integer Number of modes for which the critical bandwidth is calculated, default is 2.
#' @param BPPARAM A BiocParallelParam object specifying whether perform the analysis parallelly using 
#' \code{BiocParallel} default is \code{SerialParam()}, meaning no parallel. 
#' You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=T)} to parallel it, 
#' the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}. default is \code{SerialParam()}.
#' @param name character the assay name in \code{assay} of the result, default is \code{ModeFlag}.
#' @param features Vector specifying the subset of features to be used for the analysis, if \code{gsvaexp} is 
#' specified, it should be the rownames or rownames index of \code{gsvaexp(data)}, default is NULL, meaning
#' all features will be analyzed.
#' @param cells Vector specifying the subset of cells to be used for the analysis, default is NULL, meaning
#' all cells will be analyzed.
#' @param ... dot parameters.
#' @return if input is a \linkS4class{SVPExperiment}, output will be also a \linkS4class{SVPExperiment}, and the result
#' was stored in the \code{assays} of \linkS4class{SingleCellExperiment} of the specified \code{gsvaexp}, 
#' which is a \linkS4class{SingleCellExperiment}.  If input is a \linkS4class{SingleCellExperiment} (which is extracted from
#' \linkS4class{SVPExperiment} using \code{gsvaExp()} funtion), output will be a \linkS4class{SingleCellExperiment},
#' the result can be extracted using \code{assay()} function with specified name in default is \code{ModeFlag}.
#' @seealso to calculate the activity score of gene sets or pathway: [`runSGSA`],
#' to keep the max gene set or pathway activity score of cell: [`cluster.assign`].
#' @author Shuangbin Xu
#' @examples
#' library(SpatialExperiment)
#' # This result can be extract from the
#' # result of runSGSA with gsvaExp(svpe)
#' data(hpda_spe_cell_dec)
#' assay
#' hpda_spe_cell_dec <- hpda_spe_cell_dec |> 
#'    pred.feature.mode(assay.type = 'rwr.score', 
#'      mod0 = 2, 
#'      BPPARAM=BiocParallel::MulticoreParam(workers=2,progressbar=TRUE)
#'    )
#' assays(hpda_spe_cell_dec)
#' # We extract the activity mode of Cancer clone A and Cancer clone B, then visualize them
#' assay(hpda_spe_cell_dec, 'ModeFlag') |> t() |> 
#'  magrittr::extract(,c(2, 3)) |> 
#'  as.matrix() |> 
#'  data.frame() |> 
#'  dplyr::mutate_all(as.factor) |> 
#'  DataFrame() -> colData(hpda_spe_cell_dec)
#' \dontrun{
#'   library(ggplot2)
#'   library(ggsc)
#'   p1 <- sc_spatial(hpda_spe_cell_dec, 
#'              features = rownames(hpda_spe_cell_dec), 
#'              mapping = aes(x=x,y=y, color=Cancer.clone.A), 
#'              plot.pie = T, 
#'              pie.radius.scale = .8, 
#'              bg_circle_radius = 1.1, 
#'              color=NA, 
#'              linewidth=2
#'   ) + 
#'   scale_color_manual(values=c('white', 'black'))
#'   p1
#'   f1 <- sc_spatial(hpda_spe_cell_dec, features="Cancer clone A", 
#'              mapping=aes(x=x,y=y),
#'              pointsize=10
#'   ) + 
#'   geom_scattermore2(
#'     mapping = aes(bg_color=Cancer.clone.A, subset=Cancer.clone.A==1), 
#'     bg_line_width = .15, 
#'     pointsize = 8
#'   ) + 
#'   scale_bg_color_manual(values=c('black')) 
#'   f1
#'   
#'   f2 <- sc_spatial(hpda_spe_cell_dec, features="Cancer clone B",
#'              mapping=aes(x=x,y=y),
#'              pointsize=10
#'   ) +
#'   geom_scattermore2(
#'     mapping = aes(bg_color=Cancer.clone.B, subset=Cancer.clone.B==1),
#'     bg_line_width = .15,
#'     pointsize = 8
#'   ) +
#'   scale_bg_color_manual(values=c('black')) 
#'   f2  
#' }
setGeneric("pred.feature.mode",
  function(
    data,
    assay.type = 'rwr.score',
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    mod0 = 2,
    BPPARAM = SerialParam(),
    name = 'ModeFlag',
    features = NULL,
    cells = NULL,
    ...
  )
  standardGeneric('pred.feature.mode')
)

#' @rdname pred.feature.mode-method
#' @aliases pred.feature.mode,SingleCellExperiment
#' @export pred.feature.mode
setMethod(
    'pred.feature.mode',
    'SingleCellExperiment',
    function(
        data,
        assay.type = 'affi.score',
        gsvaexp = NULL,
        gsvaexp.assay.type = NULL,
        mod0 = 2,
        BPPARAM = SerialParam(),
        name = "ModeFlag",
        features = NULL,
        cells = NULL,
        ...
    ){

  if (mod0 < 2){
      rlang::abort("The {.var mod0} should be an integer larger than 1, default is 2")
  }
  if (is.null(assay.type)){
      assay.type <- 1
  }
  
  if (!is.null(features)){
      data <- data[features,,drop=FALSE]
  }
  if (!is.null(cells)){
      data <- data[,cells, drop=FALSE]
  }

  x <- assay(data, assay.type)

  assay(data, name) <- .internal.predict.feature.mode(x, mod0, BPPARAM, ...)

  return(data)
})


#' @rdname pred.feature.mode-method
#' @aliases pred.feature.mode,SVPExperiment
#' @export pred.feature.mode
setMethod(
    'pred.feature.mode',
    'SVPExperiment',
    function(
        data,
        assay.type = 'affi.score',
        gsvaexp = NULL,
        gsvaexp.assay.type = NULL,
        mod0 = 2,
        BPPARAM = SerialParam(),
	name = "ModeFlag",
        features = NULL,
        cells = NULL,
        ...
    ){
    
    if (!is.null(cells)){
        data <- data[,cells, drop=FALSE]
    }

    if (!is.null(gsvaexp)){
       cli::cli_inform("The {.var gsvaexp} was specified, the specified {.var gsvaExp} will be used to predict the features mode.")
 
       da2 <- gsvaExp(data, gsvaexp, withColData=FALSE, withSpatialCoords=FALSE, withImageData = FALSE)
       da2 <- pred.feature.mode(da2, assay.type = gsvaexp.assay.type, mod0 = mod0, 
                                BPPARAM = BPPARAM, name = name, features = features, cells = cells, ...)
       gsvaExp(data, gsvaexp) <- da2
    }else{
       data <- callNextMethod()
    }
    return(data)
})


.internal.predict.feature.mode <- function(da, mod0 = 2, BPPARAM = SerialParam(), ...){
    res <- BiocParallel::bplapply(seq(nrow(da)), function(x){
         .internal.locmodes(da[x,], mod0 = mod0, ...)
      },
      BPPARAM = BPPARAM
    )

    res <- do.call('rbind', res) |> Matrix::Matrix(sparse=TRUE)

    rownames(res) <- rownames(da)
    colnames(res) <- colnames(da)
    
    return(res)
}

.internal.locmodes <- function(x, mod0 = 2, ...){
    rlang::check_installed("multimode", "`pred.feature.mode()`")
    
    res <- suppressWarnings(multimode::locmodes(x,
                        mod0 = mod0,
                        display = FALSE,
                        ...))

    nm <- length(res$location) 
    if (nm %% 2 != 0){
        nm <- nm - 1 
    }
     
    index <- (seq(nm) %% 2) == 0
    
    mode.values <- res$location[index]
    res <- findIntervalCpp(x, mode.values) |> as.vector()
    return(res)

}
