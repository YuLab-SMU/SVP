#' @title One hot encode for the specified cell category.
#' @description
#' This function convert the specified cell category to one hot encode
#' @rdname runENCODE-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param group.by character a specified category column names (for example the cluster column name) of
#' \code{colData(data)}, required.
#' @param rm.group.nm character which want to remove some group type names from the names of 
#' the specified category group, default is NULL.
#' @param ... currently meaningless.
#' @examples
#' data(sceSubPbmc)
#' sceSubPbmc
#' sceSubPbmc <- runENCODE(sceSubPbmc, group.by = 'seurat_annotations') 
#' sceSubPbmc
#' gsvaExp(sceSubPbmc, 'seurat_annotations')
#' sceSubPbmc <- runENCODE(sceSubPbmc, group.by = 'seurat_annotations', rm.group.nm = c('Platelet'))
#' sceSubPbmc
#' gsvaExp(sceSubPbmc, 'seurat_annotations')
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
    rm.group.nm=NULL){

  dt <- colData(data)[[group.by]] |> as.character()
  keep.nm <- unique(dt)
  if (!is.null(rm.group.nm)){
    keep.nm <- setdiff(keep.nm, rm.group.nm)
  }
  if (length(keep.nm)==0){
     cli::cli_wran("The {.var group.by} column does not contain any elements.")
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

  data <- .sce_to_svpe(data)
  gsvaExp(data, group.by) <- SingleCellExperiment(assays=list(counts = mt))
  return(data)
})



