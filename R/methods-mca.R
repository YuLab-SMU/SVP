#' @useDynLib SVP
#' @importFrom Rcpp sourceCpp
#' @export
setGeneric('runMCA', function(data, 
			      assay.type = 'logcounts', 
			      reduction.name = "MCA", 
			      ncomponents = 50, 
			      subset.row = NULL, 
			      consider.coord = TRUE,
			      ...)
  standardGeneric('runMCA')
)

#' @export
setMethod('runMCA', 'SingleCellExperiment', 
	  function(data, 
		   assay.type = 'logcounts', 
		   reduction.name = 'MCA', 
		   ncomponents = 50, 
		   subset.row = NULL, 
		   consider.coord = TRUE,
		   ...){
  if (assay.type == 'logcounts'){
      if (!assay.type %in% assayNames(data)){
          cli::cli_abort("the 'logcounts' is not present in the assays of {.cls {class(data)}}.")
      }
  }
  
  if (!is.null(subset.row)){
      data <- data[subset.row,]
  }
  
  x <- assay(data, assay.type)

  flag.coords <- .check_element_obj(data, key = 'spatialCoords', basefun = int_colData, namefun = names)
  if (flag.coords && consider.coord){
      specoords <- .extract_element_object(data, key = 'spatialCoords', basefun = int_colData, namefun = names)
      #specoords <- log(specoords)
      #specoords[is.infinite(specoords)] <- 0
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



