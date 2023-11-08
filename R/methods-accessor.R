#' Some accessor funtions to get the internal slots of SVPExperiment
#' @name SVP-accessors 
#' @docType methods
#' @param x a \linkS4class{SVPExperiment} class.
#' @param object a \linkS4class{SVPExperiment} class.
#' @importFrom SpatialExperiment spatialCoords
#' @aliases spatialCoords,SVPExperiment-method
#' @return matrix or character or print the information of object.
#' @export
setMethod('spatialCoords', 'SVPExperiment', function(x){
    flag <- .check_element_obj(x, key = 'spatialCoords', basefun = int_colData, namefun = names)
    if (flag){
        x <- .extract_element_object(x, key = 'spatialCoords', basefun = int_colData, namefun = names)
    }else{
        x <- NULL
    }
    return(x)
})

#' @rdname SVP-accessors
#' @importFrom SpatialExperiment spatialCoordsNames
#' @aliases spatialCoordsNames,SVPExperiment-method
#' @export
setMethod('spatialCoordsNames', 'SVPExperiment', function(x){
    colnames(spatialCoords(x))
})

#' @rdname SVP-accessors
#' @importFrom SpatialExperiment imgData
#' @aliases imgData,SVPExperiment-method
#' @export
setMethod('imgData', 'SVPExperiment', function(x){
    flag <- .check_element_obj(x, key = 'imgData', basefun = int_metadata, namefun = names)
    if (flag){
        x <- .extract_element_object(x, key = 'imgData', basefun = int_metadata, namefun = names)
    }else{
        x <- NULL
    }
    return(x)    
})
