#' Run Multiple Correspondence Analysis
#' @rdname runMCA-method
#' @description Perform a Multiple Correspondence Analysis (MCA) on cells, based on
#' the expression data in a SingleCellExperiment object. It is modified based on the
#' \code{RunMCA} of \code{CelliD} with the source codes of C++.
#' @param data a SingleCellExperiment object 
#' @param assay.type which expressed data to be pulled to run, 
#' default is \code{logcounts}.
#' @param reduction.name name of the reduction result, default is \code{MCA}.
#' @param ncomponents number of components to compute and store, default is 30.
#' @param subset.row Vector specifying the subset of features to be used for 
#' dimensionality reduction. This can be a character vector of row names, 
#' an integer vector of row indices or a logical vector, default is NULL, meaning
#' all features to be used for dimensionality reduction.
#' @param subset.col Vector specifying the subset of cells to be used for
#' dimensionality reduction. This can be a character vector of column names,
#' an integer vector of column indices or a logical vector, default is NULL, meaning
#' all cells to be used for dimensionality reduction.
#' @param group.by.vars character the name(s) of covariates that harmony will remove its
#' effect on the data, default is NULL.
#' @param consider.spcoord whether consider the spatial coords as the features of data 
#' to run MCA, default is FALSE (TRUE is experimental).
#' @param ... additional parameters, see also \code{RunHarmony}.
#' @return a \linkS4class{SingleCellExperiment} and the reduction result of \code{MCA}
#' can be extracted using \code{reducedDim()} function.
#' @export
#' @examples
#' library(scuttle)
#' library(SingleCellExperiment)
#' small.sce <- mockSCE()
#' small.sce <- logNormCounts(small.sce)
#' # To improve computational efficiency, you can use RhpcBLASctl to control the number 
#' # of threads on BLAS. From example
#' # RhpcBLASctl::blas_set_num_threads(threads = 48)
#' small.sce <- runMCA(small.sce, assay.type = 'logcounts',
#'                     reduction.name = 'MCA', ncomponents = 20) 
#' # The MCA result can be extracted using reducedDim of SingleCellExperiment
#' mca.res <- reducedDim(small.sce, 'MCA')
#' mca.res |> str()
setGeneric('runMCA', function(data, 
                              assay.type = 'logcounts', 
                              reduction.name = "MCA", 
                              ncomponents = 30, 
                              subset.row = NULL, 
                              subset.col = NULL,
                              group.by.vars = NULL,
                              consider.spcoord = FALSE,
                              ...)
  standardGeneric('runMCA')
)

#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay assayNames
#' @rdname runMCA-method
#' @aliases runMCA,SingleCellExperiment
#' @export runMCA
setMethod('runMCA', 'SingleCellExperiment', 
          function(data, 
                   assay.type = 'logcounts', 
                   reduction.name = 'MCA', 
                   ncomponents = 50, 
                   subset.row = NULL,
                   subset.col = NULL, 
                   group.by.vars = NULL,
                   consider.spcoord = FALSE,
                   ...){
  if (is.numeric(assay.type)){
      assay.type <- assayNames(data)[assay.type]
  }
  if (!assay.type %in% assayNames(data)){
      cli::cli_abort("the {.var assay.type} = {assay.type} is not present in the assays of {.cls {class(data)}}.")
  }

  x <- assay(data, assay.type)

  if (!is.null(subset.row)){
      x <- x[subset.row,]
  }

  if (!is.null(subset.col)){
      x <- x[,subset.col]
  }  

  flag.coords <- .check_element_obj(data, key = 'spatialCoords', basefun = int_colData, namefun = names)
  if (flag.coords && consider.spcoord){
      if (!is.null(subset.col)){
          specoords <- .extract_element_object(data[,subset.col], key = 'spatialCoords', basefun = int_colData, namefun = names)
      }else{
          specoords <- .extract_element_object(data, key = 'spatialCoords', basefun = int_colData, namefun = names)
      }
      specoords <- .normalize.coords(specoords)
      x <- rbind(x, t(specoords))
  }
  metadata <- NULL
  if (!is.null(group.by.vars)){
      if (!is.null(subset.col)){
          metadata <- colData(data[,subset.col])[,group.by.vars,drop=FALSE] |> as.data.frame(check.names=FALSE)
      }else{
          metadata <- colData(data)[,group.by.vars,drop=FALSE] |> as.data.frame(check.names=FALSE)
      }
  }
  res.mca <- .runMCA.internal(x,
                              metadata, 
                              reduction.name = reduction.name, 
                              ncomponents = ncomponents,
                              ...
             )

  reducedDim(data, reduction.name) <- res.mca

  return(data)
})



