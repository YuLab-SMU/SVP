#' Some accessor funtions to get the internal slots of SVPExperiment
#' @name SVP-accessors
#' @param x a \linkS4class{SVPExperiment} class.
#' @param object a \linkS4class{SVPExperiment} class.
#' @param value matrix for \code{spatialCoords(object) <- value}
#' character for \code{spatialCoordsNames(object) <- value}.
#' @importFrom SpatialExperiment spatialCoords
#' @aliases 
#' spatialCoords,SVPExperiment-method
#' spatialCoordsNames,SVPExperiment-method
#' imgData,SVPExperiment-method
#' show,SVPExperiment-method
#' spatialCoordsNames<-,SVPExperiment,character-method
#' imgData<-,SVPExperiment,DataFrame-method
#' imgData<-,SVPExperiment,NULL-method
#' @return matrix or character or print the information of object
#' or a \linkS4class{SVPExperiment} object.
#' @examples
#' library(SpatialExperiment) |> suppressPackageStartupMessages()
#' library(DropletUtils) |> suppressPackageStartupMessages()
#' example(read10xVisium, echo = FALSE)
#' svpe <- as(spe, 'SVPExperiment')
#' svpe
#' spatialCoords(svpe) |> head()
NULL

#' @rdname SVP-accessors
#' @exportMethod spatialCoords
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
#' @exportMethod spatialCoordsNames
setMethod('spatialCoordsNames', 'SVPExperiment', function(x){
    colnames(spatialCoords(x))
})

#' @rdname SVP-accessors
#' @importFrom SpatialExperiment imgData
#' @exportMethod imgData
setMethod('imgData', 'SVPExperiment', function(x){
    flag <- .check_element_obj(x, key = 'imgData', basefun = int_metadata, namefun = names)
    if (flag){
        x <- .extract_element_object(x, key = 'imgData', basefun = int_metadata, namefun = names)
    }else{
        x <- NULL
    }
    return(x)    
})

#' @rdname SVP-accessors
#' @importFrom S4Vectors isEmpty
#' @exportMethod imgData<-
setReplaceMethod('imgData', c('SVPExperiment', "DataFrame"), function(x, value){
    flag <- .check_element_obj(x, key='imgData', basefun = int_metadata, namefun = names)
    if (flag){
        if (!isEmpty(value)){
            msg <- .imgData_validity(value)
            if (!is.null(msg)){
                cli::cli_abort(msg)
            }
        }
        int_metadata(x)$imgData <- value
    }
    return(x)
})

#' @rdname SVP-accessors
#' @exportMethod imgData<-
setReplaceMethod('imgData', c("SVPExperiment", "NULL"), function(x, value){
    flag <- .check_element_obj(x, key='imgData', basefun = int_metadata, namefun = names)
    if (flag){
        value <- DataFrame()
        `imgData<-`(x, value)
    }
    return(x)
})

#' @rdname SVP-accessors
#' @aliases spatialCoords<-,SVPExperiment
#' @exportMethod spatialCoords<-
setReplaceMethod("spatialCoords", c("SVPExperiment", "matrix_Or_NULL"), function(x, value){
    flag <- .check_element_obj(x, key='spatialCoords', basefun = int_colData, namefun = names)
    if (!flag || is.null(value)){
        int_colData(x)$spatialCoords <- matrix(numeric(), ncol(x),0)
        if (is.null(value)){
            return(x)
        }
    }
    flag1 <- is.numeric(value)
    flag2 <- nrow(value) >= ncol(x)
    flag3 <- all(colnames(x) %in% rownames(value))
    if (!flag1){
       cli::cli_abort("The `value` (coordinate of cell or spot) must be a numeric matrix.")
    }
    if (!flag2){
       cli::cli_abort("The row number of coordinate matrix must be the same to the column number of {.cls class(x)}.")
    }
    if (!flag3){
       cli::cli_abort("The rownames of coordinate matrix must be the same of the column names of {.cl class(x)}.")
    }
    int_colData(x)$spatialCoords <- value[colnames(x),,drop=FALSE]
    return(x)
})


#' @rdname SVP-accessors
#' @importFrom SpatialExperiment spatialCoordsNames<-
#' @exportMethod spatialCoordsNames<-
setReplaceMethod("spatialCoordsNames", c('SVPExperiment', 'character'),
    function(x, value){
    flag <- .check_element_obj(x, key='spatialCoords', basefun = int_colData, namefun = names)
    if (flag){
        colnames(int_colData(x)$spatialCoords) <- value
    }    
    return(x)
})


.imgData_validity <- utils::getFromNamespace(".imgData_validity", "SpatialExperiment")
