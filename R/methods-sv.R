#' @title spatial or single cell variable features matrix extract method
#' 
#' @description
#' the identification result of the spatial variable or single cell variable (SV) features 
#' is important to the downstream analysis.
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
#' \describe{
#' \item{\code{svDf(x, type)}:}{
#' Retrieves a \linkS4class{DataFrame} containing the new features (gene sets) (rows) 
#' for the specified \code{type}.
#' \code{type} should either be a string specifying the name of the features scores matrix 
#' in \code{x} to retrieve, or a numeric scalar specifying the index of the desired matrix, 
#' defaulting to the first matrix is missing.
#' }
#'
#' \item{\code{svDfNames(x)}:}{
#' Retures a character vector containing the names of all features SV DataFrame Lists in
#' \code{x}. This is guaranteed to be of the same length as the number of results.
#' }
#' 
#' \item{\code{svDfs(x)}:}{
#' Returns a named \linkS4class{List} of matrices containing one or more \linkS4class{DataFrame} objects.
#' Each object is guaranteed to have the same number of rows, in a 1:1 correspondence to those in \code{x}.
#' }
#'}
#'
#' @section Single-object setter:
#' \code{svDf(x, type) <- value} will add or replace an SV matrix in a 
#' \linkS4class{SingleCellExperiment} object \code{x}.
#' The value of \code{type} determines how the result is added or replaced:
#' \itemize{
#' \item If \code{type} is missing, \code{value} is assigned to the first result.
#' If the result already exists, its name is preserved; otherwise it is given a default name \code{"unnamed.sv1"}.
#' \item If \code{type} is a numeric scalar, it must be within the range of existing results, and \code{value} will 
#' be assigned to the result at that index.
#' \item If \code{type} is a string and a result exists with this name, \code{value} is assigned to to that result.
#' Otherwise a new result with this name is append to the existing list of results.
#' }
#' 
#' @section Other setter:
#' \describe{
#' \item{\code{svDfs(x) <- value}:}{
#' Replaces all features sv result matrixs in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of \linkS4class{DataFrame} objects
#' with number of row equal to \code{nrow(x)}.
#'
#' If \code{value} is named, those names will be used to name the SV matrixs in \code{x}.
#' Otherwise, unnamed results are assigned default names prefixed with \code{"unnamed.sv"}.
#'
#' If \code{value} is \code{NULL}, all SV matrixs in \code{x} are removed.
#' }
#'
#' \item{\code{svDfNames(x) <- value}:}{
#' Replaces all names for SV matrixs in \code{x} with a character vector \code{value}.
#' This should be of length equal to the number of results currently in \code{x}.
#'}
#'}
#' @name svDfs
#' @docType methods
#' @aliases
#' svDf svDfs svDfNames
#' svDf,SingleCellExperiment,missing-method
#' svDf,SingleCellExperiment,numeric-method
#' svDf,SingleCellExperiment,character-method
#' svDfs,SingleCellExperiment-method
#' svDfNames,SingleCellExperiment-method
#' svDf<- svDfs<- svDfNames<-
#' svDf<-,SingleCellExperiment,missing-method
#' svDf<-,SingleCellExperiment,numeric-method
#' svDf<-,SingleCellExperiment,character-method
#' svDfs<-,SingleCellExperiment-method
#' svDfNames<-,SingleCellExperiment,character-method
NULL

#' @export
setGeneric('svDfs', function(x,...)standardGeneric('svDfs'))

#' @export
setGeneric('svDf', function(x, type, ...)standardGeneric('svDf'))

#' @export
setGeneric('svDfNames', function(x, ...)standardGeneric('svDfNames'))

#' @export
setGeneric("svDf<-", function(x, type, ..., value) standardGeneric("svDf<-"))

#' @export
setGeneric("svDfs<-", function(x, ..., value) standardGeneric("svDfs<-"))


#' @export
setGeneric('svDfNames<-', function(x, value)standardGeneric('svDfNames<-'))

#' @export
setMethod('svDfNames', 'SingleCellExperiment', function(x){
    colnames(int_elementMetadata(x)[[.sv_key]])
}) 


#' @export
setReplaceMethod('svDfNames', c('SingleCellExperiment', 'character'), function(x, value){
    .set_internal_names(x, value,
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.sv_key,
        unname.key = .unnamed.sv
    )
})

#' @export
setMethod("svDf", c("SingleCellExperiment", "missing"), function(x, type) {
    .get_internal_missing(x,
        basefun=svDf,
        namefun=svDfNames,
        funstr="svDf"
    )
})

#' @export
setMethod("svDf", c("SingleCellExperiment", "numeric"), function(x, type) {
    out <- .get_internal_numeric(x, type,
        getfun=int_elementMetadata,
        key=.sv_key,
        funstr="svDf",
        substr="type")

    out
})

#' @export
setMethod("svDf", c("SingleCellExperiment", "character"), function(x, type) {
    out <- .get_internal_character(x, type,
        getfun=int_elementMetadata,
        key=.sv_key,
        funstr="svDf",
        substr="type",
        namestr="svNames")

    out
})


#' @export
setReplaceMethod("svDf", c("SingleCellExperiment", "missing"), function(x, type, ..., value) {
    .set_internal_missing(x, value,
        basefun=`svDf<-`,
        namefun=svDfNames
    )
})


#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
#' @export
setReplaceMethod("svDf", c("SingleCellExperiment", "numeric"), function(x, type, ..., value) {
    .set_internal_numeric(x, type, value,
        getfun = int_elementMetadata,
        setfun = `int_elementMetadata<-`,
        key = .sv_key,
	convertfun = NULL,
	xdimfun = nrow,
	vdimfun = nrow,
        funstr = "svDf",
	xdimstr = 'nrow',
	vdimstr = 'nrow',
        substr = "type")
})


#' @export
setReplaceMethod("svDf", c("SingleCellExperiment", "character"), function(x, type, ..., value) {
    .set_internal_character(x, type, value,
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.sv_key,
        convertfun=NULL,
        xdimfun=nrow,
        vdimfun=nrow,
        funstr="svDf",
        xdimstr="nrow",
        vdimstr="nrow",
        substr="type")
})

#' @export
setMethod("svDfs", "SingleCellExperiment", function(x) {
    value <- .get_internal_all(x, 
        getfun=int_elementMetadata, 
        key=.sv_key)
    value
})


#' @export
setReplaceMethod("svDfs", "SingleCellExperiment", function(x, value) {
    .set_internal_all(x, value,
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.sv_key,
        convertfun= NULL,
        xdimfun=nrow,
        vdimfun=nrow,
        funstr="svDfs",
        xdimstr="nrow",
        vdimstr="nrow")
})
