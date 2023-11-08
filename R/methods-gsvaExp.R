#' Gene Set Variation Analysis Experiment methods
#'
#' @description
#' In some experiments, gene set variation analysis will generated different features (the names of KEGG pathway or the GO term).
#' These data cannot be stored in the main \code{assays} of the \linkS4class{SVPExperiment} itself.
#' However, it is still desirable to store these features \emph{somewhere} in the SVPExperiment.
#' This simplifies book-keeping in long workflows and ensure that samples remain synchronised.
#'
#' To facilitate this, the \linkS4class{SVPExperiment} class allows for \dQuote{gene set variation analysis experiments}.
#' Nested \linkS4class{SingleCellExperiment}-class objects are stored inside the SVPExperiment object \code{x}, 
#' in a manner that guarantees that the nested objects have the same columns in the same order as those in \code{x}.
#' Methods are provided to enable convenient access to and manipulation of these gene set variation analysis Experiments.
#' Each GSVA Experiment should contain experimental data and row metadata for a distinct set of features.
#' (These methods refer to the \code{altExp} of \code{SingleCellExperiment}).
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SVPExperiment} object.
#' \describe{
#' \item{\code{gsvaExp(x, e, withDimnames=TRUE, withColData=TRUE, withSpatialCoords = TRUE, withImgData=TRUE, withReducedDim=FALSE)}:}{
#' Retrieves a \linkS4class{SingleCellExperiment} containing gene set name features (rows) for all cells (columns) in \code{x}.
#' \code{e} should either be a string specifying the name of the gene set variation Experiment in \code{x} to retrieve,
#' or a numeric scalar specifying the index of the desired Experiment, defaulting to the first Experiment is missing.
#'
#' \code{withDimnames=TRUE}, the column names of the output object are set to \code{colnames(x)}.
#' In addition, if \code{withColData=TRUE}, \code{\link{colData}(x)} is \code{cbind}ed to the front of the column data of the output object.
#' \code{withSpatialCoords = TRUE}, the spatial coordinates of the output object are set to \code{spatialCoords(x)} if \code{x} has 
#' spatial coordinates.
#' \code{withImgData=TRUE}, the image metadata of the output object are set to \code{imgData(x)} if \code{x} has image metadata.
#' If \code{withReducedDim=TRUE}, the dimensionality reduction results of output object are set to \code{reducedDims(x)} if \code{x} has 
#' dimensionality reduction results
#' }
#' 
#' \item{\code{gsvaExpNames(x)}:}{
#' Returns a character vector containing the names of all gene set variation Experiments in \code{x}.
#' This is guaranteed to be of the same length as the number of results, though the names may not be unique.
#' }
#' \item{\code{gsvaExps(x, withDimnames=TRUE, withColData=TRUE, withSpatialCoords = TRUE, withImgData=TRUE, withReducedDim=FALSE)}:}{
#' Returns a named \linkS4class{List} of matrices containing one or more \linkS4class{SingleCellExperiment} objects.
#' Each object is guaranteed to have the same number of columns, in a 1:1 correspondence to those in \code{x}.
#'
#' If \code{withDimnames=TRUE}, the column names of each output object are set to \code{colnames(x)}.
#' In addition, if \code{withColData=TRUE}, \code{\link{colData}(x)} is \code{cbind}ed to the front of the column data of each output object.
#' \code{withSpatialCoords = TRUE}, the spatial coordinates of the output object are set to \code{spatialCoords(x)} if \code{x} has
#' spatial coordinates.
#' \code{withImgData=TRUE}, the image metadata of the output object are set to \code{imgData(x)} if \code{x} has image metadata.
#' If \code{withReducedDim=TRUE}, the dimensionality reduction results of output object are set to \code{reducedDims(x)} if \code{x} has
#' dimensionality reduction results
#' }
#' }
#'
#' @section Single-object setter:
#' \code{gsvaExp(x, e, withDimnames=TRUE, withColData=FALSE, withSpatialCoords = FALSE, withImgData = FALSE, withReducedDim = FALSE) <- value} will 
#' add or replace an gene set variation Experiment in a \linkS4class{SVPExperiment} object \code{x}.
#' The value of \code{e} determines how the result is added or replaced:
#' \itemize{
#' \item If \code{e} is missing, \code{value} is assigned to the first result.
#' If the result already exists, its name is preserved; otherwise it is given a default name \code{"unnamed.gsva1"}.
#' \item If \code{e} is a numeric scalar, it must be within the range of existing results, and \code{value} will be assigned to the result at that index.
#' \item If \code{e} is a string and a result exists with this name, \code{value} is assigned to to that result.
#' Otherwise a new result with this name is append to the existing list of results.
#' }
#'
#' \code{value} is expected to be a SingleCellExperiment object with number of columns equal to \code{ncol(x)}.
#' Alternatively, if \code{value} is \code{NULL}, the gene set variation Experiment at \code{e} is removed from the object.
#'
#' If \code{withDimnames=TRUE}, the column names of \code{value} are checked against those of \code{x}.
#' A warning is raised if these are not identical, with the only exception being when \code{value=NULL}.
#' This is inspired by the argument of the same name in \code{\link{assay<-}}.
#'
#' If \code{withColData=TRUE}, we assume that the left-most columns of \code{colData(value)} are identical to \code{colData(x)}.
#' If so, these columns are removed, effectively reversing the \code{withColData=TRUE} setting for the \code{gsvaExp} getter.
#' Otherwise, a warning is raised.
#' 
#' If \code{withSpatialCoords = TRUE}, the spatial coordinates will be kept in the \code{value} if it has, and will add or replace it in a 
#' \linkS4class{SVPExperiment} object \code{x}.
#' 
#' If \code{withImgData = TRUE}, the image metadata will be kept in the \code{value} if it has, and will add or replace it in a 
#' \linkS4class{SVPExperiment} object \code{x}.
#'
#' If \code{withReducedDim = TRUE}, the dimensionality reduction results will be kept in the \code{value} if it has, and will add or replace it
#' in a \linkS4class{SVPExperiment} object \code{x}.
#'
#'
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SVPExperiment} object.
#' \describe{
#' \item{\code{gsvaExps(x, withDimnames=TRUE, withColData=FALSE, withSpatialCoords = FALSE, withImgData = FALSE, withReducedDim = FALSE) <- value}:}{
#' Replaces all gene set variation Experiments in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of SingleCellExperiment objects
#' with number of columns equal to \code{ncol(x)}.
#'
#' If \code{value} is named, those names will be used to name the gene set variant Experiments in \code{x}.
#' Otherwise, unnamed results are assigned default names prefixed with \code{"unnamed.gsva"}.
#'
#' If \code{value} is \code{NULL}, all gene set variation Experiments in \code{x} are removed.
#'
#' If \code{value} is a \linkS4class{Annotated} object, any \code{\link{metadata}} will be retained in \code{gsvaExps(x)}.
#' If \code{value} is a \linkS4class{Vector} object, any \code{\link{mcols}} will also be retained.
#'
#' If \code{withDimnames=TRUE}, the column names of each entry of \code{value} are checked against those of \code{x}.
#' A warning is raised if these are not identical.
#'
#' If \code{withColData=TRUE}, we assume that the left-most columns of the \code{colData} for each entry of \code{value} are identical to \code{colData(x)}.
#' If so, these columns are removed, effectively reversing the \code{withColData=TRUE} setting for the \code{gsvaExps} getter.
#' Otherwise, a warning is raised.
#' If \code{withSpatialCoords = TRUE}, \code{withImgData = TRUE}, and \code{withReducedDim = TRUE} refer to the \code{gsvaExp(...) <- value}.
#' }
#' \item{\code{gsvaExpNames(x) <- value}:}{
#' Replaces all names for gene set variant Experiments in \code{x} with a character vector \code{value}.
#' This should be of length equal to the number of results currently in \code{x}.
#' }
#' }
#'
#'
#' @section Main Gene Set Variation Experiment naming:
#' The Gene Set Variation Experiments are naturally associated with names (\code{e} during assignment).
#' However, we can also name the main Experiment in a \linkS4class{SVPExperiment} \code{x}:
#' \describe{
#' \item{\code{mainGsvaExpName(x) <- value}:}{
#' Set the name of the main Experiment to a non-\code{NA} string \code{value}.
#' This can also be used to unset the name if \code{value=NULL}.
#' }
#' \item{\code{mainGsvaExpName(x)}:}{
#' Returns a string containing the name of the main Experiment.
#' This may also be \code{NULL} if no name is specified.
#' }
#' }
#' 
#'
#' @name gsvaExps
#' @docType methods
#' @aliases
#' gsvaExp gsvaExps gsvaExpNames
#' gsvaExp,SVPExperiment,missing-method
#' gsvaExp,SVPExperiment,numeric-method
#' gsvaExp,SVPExperiment,character-method
#' gsvaExps,SVPExperiment-method
#' gsvaExpNames,SVPExperiment-method
#' gsvaExp<- gsvaExps<- gsvaExpNames<-
#' gsvaExp<-,SVPExperiment,missing-method
#' gsvaExp<-,SVPExperiment,numeric-method
#' gsvaExp<-,SVPExperiment,character-method
#' gsvaExps<-,SVPExperiment-method
#' gsvaExpNames<-,SVPExperiment,character-method
#' [,SCEByColumn,ANY,ANY,ANY-method
#' [<-,SCEByColumn,ANY,ANY,ANY-method
#' c,SCEByColumn-method
#' length,SCEByColumn-method
#' names,SCEByColumn-method
#' names<-,SCEByColumn-method
#' mainGsvaExpName
#' mainGsvaExpName,SVPExperiment-method
#' mainGsvaExpName<-
#' mainGsvaExpName<-,SVPExperiment,character_OR_NULL-method
#' @return see \code{Getter} and \code{setter}.
NULL

#' @export
setGeneric('gsvaExps', function(x,...)standardGeneric('gsvaExps'))

#' @export
setGeneric('gsvaExp', function(x, e, withDimnames=TRUE, withColData=TRUE, withSpatialCoords = TRUE, withImgData=TRUE, withReducedDim=FALSE, ...)standardGeneric('gsvaExp'))

#' @export
setGeneric('gsvaExpNames', function(x, ...)standardGeneric('gsvaExpNames'))

#' @export
setGeneric('mainGsvaExpName', function(x)standardGeneric('mainGsvaExpName'))

#' @export
setGeneric("gsvaExp<-", function(x, e, withDimnames=TRUE, withColData=FALSE, withSpatialCoords = FALSE, withImgData=FALSE, withReducedDim=FALSE, ..., value) standardGeneric("gsvaExp<-"))

#' @export
setGeneric("gsvaExps<-", function(x, withDimnames=TRUE, withColData=FALSE, withSpatialCoords = FALSE, withImgData = FALSE, withReducedDim=FALSE, ..., value) standardGeneric("gsvaExps<-"))

#' @export
setGeneric('mainGsvaExpName<-', function(x, value)standardGeneric('mainGsvaExpName<-'))

#' @export
setGeneric('gsvaExpNames<-', function(x, value)standardGeneric('gsvaExpNames<-'))

#' @export
setMethod('gsvaExpNames', 'SVPExperiment', function(x){
    colnames(int_colData(x)[[.gsva_key]])
}) 


#' @export
setReplaceMethod('gsvaExpNames', c('SVPExperiment', 'character'), function(x, value){
    .set_internal_names(x, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.gsva_key
    )
})

#' @importFrom SingleCellExperiment int_metadata<- int_metadata
#' @export
setMethod("mainGsvaExpName", "SVPExperiment", function(x){
    int_metadata(x)$mainGsvaExpName
})

#' @export
setReplaceMethod("mainGsvaExpName", c("SVPExperiment", "character_OR_NULL"), function(x, value){
    int_metadata(x)$mainGsvaExpName <- value
    x
})

#' @export
setReplaceMethod('gsvaExp', c('SVPExperiment', 'missing'), function(x, e, withDimnames = TRUE, withColData = FALSE, withSpatialCoords=FALSE, withImgData=FALSE, withReducedDim=FALSE, ..., value){
    .set_internal_missing(x, value,
        withDimnames=withDimnames,
        withColData=withColData,
	withSpatialCoords = withSpatialCoords,
	withImgData = withImgData,
	withReducedDim = withReducedDim,
        basefun=`gsvaExp<-`,
        namefun=gsvaExpNames
    )

})

#' @export
setReplaceMethod('gsvaExp', c('SVPExperiment', 'character'), function(x, e, withDimnames = TRUE, withColData = FALSE, withSpatialCoords=FALSE, withImgData=FALSE, withReducedDim=FALSE, ..., value){
    value <- .check_gsvaexp_columns(x, value, withDimnames=withDimnames, withColData=withColData, 
				    withSpatialCoords = withSpatialCoords, withImgData = withImgData, 
				    withReducedDim = withReducedDim)
    .set_internal_character(x, e, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.gsva_key,
        convertfun=SCEByColumn,
        xdimfun=ncol,
        vdimfun=length,
        funstr='gsvaExp',
        xdimstr="ncol",
        vdimstr="columns",
        substr="e")

})


#' @importFrom SingleCellExperiment int_colData<-
#' @export
setReplaceMethod('gsvaExp', c('SVPExperiment', 'numeric'), function(x, e, withDimnames = TRUE, withColData = FALSE, withSpatialCoords=FALSE, withImgData=FALSE, withReducedDim=FALSE, ..., value){
    value <- .check_gsvaexp_columns(x, value, withDimnames=withDimnames, withColData=withColData, 
				    withSpatialCoords = withSpatialCoords, withImgData = withImgData, 
				    withReducedDim = withReducedDim)

    .set_internal_numeric(x, e, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.gsva_key,
        convertfun=SCEByColumn,
        xdimfun=ncol,
        vdimfun=length,
        funstr="gsvaExp",
        xdimstr="ncol",
        vdimstr="columns",
        substr="e")

})


#' @export
setReplaceMethod('gsvaExps', 'SVPExperiment', function(x, withDimnames=TRUE, withColData=FALSE, withSpatialCoords=FALSE, withImgData=FALSE, withReducedDim=FALSE, ..., value){
    if (withDimnames || withColData || withSpatialCoords || withImgData || withReducedDim) {
        for (v in seq_along(value)) {
            value[[v]] <- .check_gsvaexp_columns(x, value[[v]],
                withDimnames=withDimnames, withColData=withColData, 
		withSpatialCoords = withSpatialCoords, withImgData = withImgData,  
		withReducedDim = withReducedDim,
                fun='gsvaExps', vname=sprintf("value[[%s]]", v))
        }
    }

    .set_internal_all(x, value,
        getfun=int_colData,
        setfun=`int_colData<-`,
        key=.gsva_key,
        convertfun=SCEByColumn,
        xdimfun=ncol,
        vdimfun=length,
        funstr="gsvaExps",
        xdimstr="ncol",
        vdimstr="columns")
})

#' @importFrom SingleCellExperiment int_colData
#' @importFrom S4Vectors endoapply
#' @export
setMethod('gsvaExps', 'SVPExperiment', function(x, withDimnames=TRUE, withColData=TRUE, withSpatialCoords=TRUE, 
						withImgData=TRUE, withReducedDim=FALSE, ...){
    y <- .get_internal_all(x,
           getfun=int_colData,
           key=.gsva_key
         )

    y <- endoapply(y, .get_sce)

    if (withDimnames || withColData){
        y <- endoapply(y, .fill_gsvaexps_info, 
                       x = x, 
                       withDimnames = withDimnames, 
                       withColData = withColData,
                       withSpatialCoords = withSpatialCoords, 
                       withImgData = withImgData,
		       withReducedDim = withReducedDim
	     )
    }
    return(y)
})

#' @export
setMethod('gsvaExp', c('SVPExperiment', "missing"), function(x, e, withDimnames = TRUE, withColData = TRUE, 
							     withSpatialCoords = TRUE, withImgData = TRUE, 
							     withReducedDim=FALSE, ...){
    y <- .get_internal_missing(x,
            basefun = gsvaExp,
            namefun = gsvaExpNames,
            funstr = 'gsvaExp',
            withDimnames = withDimnames,
            withColData = withColData,
            withSpatialCoords = withSpatialCoords,
            withImgData = withImgData,
	    withReducedDim = withReducedDim,
            ...
         )
    return(y)
})

#' @export
setMethod('gsvaExp', c('SVPExperiment', 'numeric'), function(x, e, withDimnames = TRUE, withColData = TRUE, withSpatialCoords = TRUE, 
							     withImgData = TRUE, withReducedDim=FALSE, ...){
    y <- .get_internal_numeric(
           x,
           index = e,
           getfun = int_colData,
           key = .gsva_key,
           funstr = 'gsvaExp',
           substr = 'e'
         )
    y <- .get_sce(y)

    y <- .fill_gsvaexps_info(y, x, withDimnames, withColData, withSpatialCoords, withImgData, withReducedDim)
    return(y)

})


#' @export
setMethod('gsvaExp', c('SVPExperiment', 'character'), function(x, e, withDimnames = TRUE, withColData = TRUE, 
							       withSpatialCoords = TRUE, withImgData = TRUE, 
							       withReducedDim=FALSE, ...){
    y <- .get_internal_character(
           x,
           index = e,
           getfun = int_colData,
           key = .gsva_key,
           funstr = 'gsvaExp',
           substr = 'e',
           namestr = 'gsvaExpNames'
         )
    y <- .get_sce(y)
    y <- .fill_gsvaexps_info(y, x, withDimnames, withColData, withSpatialCoords, withImgData, withReducedDim)
    return(y)
})

