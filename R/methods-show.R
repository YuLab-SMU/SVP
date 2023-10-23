#' @importFrom methods as callNextMethod is new
#' @importFrom S4Vectors coolcat
.show_svpe <- function(object){
    callNextMethod()
    coolcat("spatialCoords names(%d) : %s\n", spatialCoordsNames(object))
    coolcat("imgData names(%d): %s\n", names(imgData(object))) 
    coolcat("gsvaExps names(%d) : %s\n", gsvaExpNames(object))
}

#' @rdname SVP-accessors
#' @importFrom methods show
#' @aliases show,SVPExperiment-method
#' @export
setMethod('show', 'SVPExperiment', .show_svpe)
