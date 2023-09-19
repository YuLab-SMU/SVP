#' @useDynLib SVP
#' @importFrom Rcpp sourceCpp
#' @export
setGeneric('runMCA', function(data, 
			      assay.type = 'logcounts', 
			      reduction.name = "MCA", 
			      ncomponents = 50, 
			      subset.row = NULL, 
			      ...)
  standardGeneric('runMCA')
)

#' @export
setMethod('runMCA', 'SingleCellExperiment', 
	  function(data, 
		   assay.type = 'logcounts', 
		   reduction.name = 'MCA', 
		   ncomponents = 50, 
		   subset.row = NULL, ...){
  if (assay.type == 'logcounts'){
      if (!assay.type %in% assayNames(data)){
          cli::cli_abort("the 'logcounts' is not present in the assays of {.cls {class(data)}}.")
      }
  }
  
  if (!is.null(subset.row)){
      data <- data[subset.row,]
  }
  
  x <- assay(data, assay.type)
  
  res.mca <- .runMCA.internal(x, 
			      reduction.name = reduction.name, 
			      ncomponents = ncomponents
             )

  reducedDim(data, reduction.name) <- res.mca

  return(data)
})



