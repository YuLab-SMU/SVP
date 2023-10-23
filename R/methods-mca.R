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
#' @param consider.spcoord whether consider the spatial coords to run MCA with 
#' the features of data, default is TRUE.
#' @param ... additional parameters, meaningless now.
#' @export
setGeneric('runMCA', function(data, 
			      assay.type = 'logcounts', 
			      reduction.name = "MCA", 
			      ncomponents = 30, 
			      subset.row = NULL, 
			      consider.spcoord = TRUE,
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
		   consider.spcoord = TRUE,
		   ...){
  if (!assay.type %in% assayNames(data)){
      cli::cli_abort("the {.var assay.type} = {assay.type} is not present in the assays of {.cls {class(data)}}.")
  }
  
  if (!is.null(subset.row)){
      data <- data[subset.row,]
  }
  
  x <- assay(data, assay.type)

  flag.coords <- .check_element_obj(data, key = 'spatialCoords', basefun = int_colData, namefun = names)
  if (flag.coords && consider.spcoord){
      specoords <- .extract_element_object(data, key = 'spatialCoords', basefun = int_colData, namefun = names)
      specoords <- .normalize.coords(specoords)
      x <- rbind(x, t(specoords))
  }
  
  res.mca <- .runMCA.internal(x, 
                              reduction.name = reduction.name, 
			      ncomponents = ncomponents
             )

  reducedDim(data, reduction.name) <- res.mca

  return(data)
})



