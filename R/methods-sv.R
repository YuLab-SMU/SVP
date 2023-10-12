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
