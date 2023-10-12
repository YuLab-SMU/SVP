#' @importFrom SpatialExperiment spatialCoords
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

#' @importFrom SpatialExperiment spatialCoordsNames
#' @export
setMethod('spatialCoordsNames', 'SVPExperiment', function(x){
    colnames(spatialCoords(x))
})

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
