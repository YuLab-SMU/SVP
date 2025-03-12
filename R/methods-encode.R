#' @title One hot encode for the specified cell category.
#' @description
#' This function convert the specified cell category to one hot encode
#' @rdname runENCODE-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param group.by character a specified category column names (for example the cluster column name) of
#' \code{colData(data)}. Or a vector of length equal to ‘ncol(data)’, specifying the group to which each cell
#' is assigned. It is required.
#' @param rm.group.nm character which want to remove some group type names from the names of 
#' the specified category group, default is NULL.
#' @param ... currently meaningless.
#' @return SVPExperiment object
#' @examples
#' data(sceSubPbmc)
#' sceSubPbmc
#' sceSubPbmc <- runENCODE(sceSubPbmc, group.by = 'seurat_annotations') 
#' sceSubPbmc
#' gsvaExp(sceSubPbmc, 'seurat_annotations')
#' sceSubPbmc <- runENCODE(sceSubPbmc, group.by = 'seurat_annotations', rm.group.nm = c('Platelet'))
#' sceSubPbmc
#' gsvaExp(sceSubPbmc, 'seurat_annotations')
#' # The group.by also can be a vector of length equal to ncol(data).
#' sceSubPbmc <- runENCODE(
#'                 sceSubPbmc, 
#'                 group.by = sceSubPbmc$seurat_annotations, 
#'                 rm.group.nm = c('Platelet')
#'               )
#' sceSubPbmc
#' identical(gsvaExp(sceSubPbmc, 'seurat_annotations'), gsvaExp(sceSubPbmc, "ENCODE"))
setGeneric('runENCODE',
  function(
    data,
    group.by,
    rm.group.nm = NULL,
    ...
  )
  standardGeneric('runENCODE')
)


#' @rdname runENCODE-method
#' @aliases runENCODE,SingleCellExperiment
#' @export runENCODE
setMethod("runENCODE", "SingleCellExperiment", 
  function(
    data, 
    group.by, 
    rm.group.nm=NULL,
    ...){

  mt <- .build_assays_from_group(data, group.by, rm.group.nm)
  data <- .sce_to_svpe(data)
  if (length(group.by) > 1){
    group.by <- 'ENCODE'
  }
  gsvaExp(data, group.by) <- SingleCellExperiment(assays=list(counts = mt))
  return(data)
})


#' @title calculate the F1 value based on LISA result in the specified category.
#' @description
#' this is to calculate the F1 value based on LISA result in some spatial domain.
#' If a feature has a larger F1 value in a spatial domain, it means the feature 
#' is more concentrated in that spatial domain (specified category). 
#' @rdname cal_lisa_f1-method
#' @param data a \linkS4class{SingleCellExperiment} object
#' @param lisa.res list the result of \code{runLISA} or \code{runLOCALBV}.
#' @param group.by character a specified category column names (for example the cluster column name) of
#' \code{colData(data)}. Or a vector of length equal to \code{ncol(data)}, specifying the group to 
#' which each cell is assigned. It is required.
#' @param type character the type of \code{cluster.test} column of result of \code{runLISA} or
#' \code{runLOCALBV}, default is \code{'High'}.
#' @param rm.group.nm character which want to remove some group type names from the names of
#' the specified category group, default is NULL.
#' @param ... currently meaningless.
#' @return a data.frame object containing the F1 value for each category in \code{group.by}.
#' @examples
#' data(hpda_spe_cell_dec)
#' lisa.res1 <- hpda_spe_cell_dec |> 
#'    runLISA(
#'      features = rownames(hpda_spe_cell_dec), 
#'      assay.type = 1
#'    )
#' res <- cal_lisa_f1(hpda_spe_cell_dec, lisa.res1, type='High', group.by = 'cluster_domain')
#' head(res)
#' # group.by, a vector of length equal to the ncol(data).
#' res2 <- cal_lisa_f1(hpda_spe_cell_dec, 
#'                     lisa.res1, 
#'                     type='High', 
#'                     group.by = hpda_spe_cell_dec$cluster_domain
#'   )
#' identical(res, res2)
setGeneric("cal_lisa_f1", 
  function(
    data,
    lisa.res,
    type = 'High',
    group.by,
    rm.group.nm = NULL,
    ...
  )
  standardGeneric("cal_lisa_f1")
)

#' @rdname cal_lisa_f1-method
#' @aliases cal_lisa_f1,SingleCellExperiment
#' @export cal_lisa_f1
setMethod("cal_lisa_f1", "SingleCellExperiment",
  function(
  data,
  lisa.res,
  type = 'High',
  group.by, 
  rm.group.nm = NULL, 
  ...){

  if (missing(lisa.res)){
    if (is.null(data@int_colData$localResults)){
       cli::cli_warn(c("The {.var lisa.res} should be provided, when the {.var data} do not contain ",
                       "the result of runLISA or runLOCALBV (speficied action='add')."))
    }
    lisa.res <- data@int_colData$localResults |> lapply(function(x) x|> as.list() |> SimpleList())
    lisa.res <- lisa.res[[1]]  
  }
  
  mt <- .build_assays_from_group(data, group.by, rm.group.nm)
  x <- lapply(lisa.res, function(x)x[,'cluster.test',drop=FALSE]) 
  x <- do.call('cbind', x) |> setNames(names(x)) |> as.matrix()
  flag <- x == type
  x[flag] <- 1
  x[!flag] <- 0
  class(x) <- 'numeric'
  x <- Matrix::Matrix(as.matrix(x), sparse=TRUE)  
  if (!all(rownames(x)==colnames(mt))){
    cli::cli_abort(c("The column names of {.var data} should be equal the rownames of ",
                     "each data frame in {.var lisa.res}."))
  }
  res <- CalF1Parallel(x, mt)
  rownames(res) <- colnames(x)
  colnames(res) <- rownames(mt)
  return(res)
})


.build_assays_from_group <- function(data, group.by, rm.group.nm){
  if (.check_group.by(group.by, ncol(data))){
    dt <- group.by |> as.character()
  }else{
    dt <- colData(data)[[group.by]] |> as.character()
  }
  keep.nm <- unique(dt)
  if (!is.null(rm.group.nm)){
    keep.nm <- setdiff(keep.nm, rm.group.nm)
  }
  if (length(keep.nm)==0){
     cli::cli_warn("The {.var group.by} column does not contain any elements.")
     return(data)
  }

  ind.j <- lapply(keep.nm, function(x) which(dt==x))
  ind.i <- lapply(seq(length(ind.j)),function(i)rep(i,each=length(ind.j[[i]]))) |> unlist()
  mt <- Matrix::sparseMatrix(
          i= ind.i,
          j = ind.j |> unlist(),
          x = 1,
          dims = c(length(keep.nm), ncol(data))
  )
  colnames(mt)<- colnames(data)
  rownames(mt) <- keep.nm
  return(mt)
}


