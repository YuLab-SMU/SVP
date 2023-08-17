#' @importFrom S4Vectors coolcat
.show_svpe <- function(object){
    callNextMethod()
    coolcat("spatialCoords names(%d) : %s\n", spatialCoordsNames(object))
    coolcat("imgData names(%d): %s\n", names(imgData(object))) 
    coolcat("gsvaExps names(%d) : %s\n", gsvaExpNames(object))
}

#' @export
setMethod('show', 'SVPExperiment', .show_svpe)
