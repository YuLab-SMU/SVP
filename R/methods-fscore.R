#' @export
setGeneric('fscoreDfs', function(x,...)standardGeneric('fscoreDfs'))

#' @export
setGeneric('fscoreDf', function(x, type, ...)standardGeneric('fscoreDf'))

#' @export
setGeneric('fscoreDfNames', function(x, ...)standardGeneric('fscoreDfNames'))

#' @export
setGeneric("fscoreDf<-", function(x, type, ..., value) standardGeneric("fscoreDf<-"))

#' @export
setGeneric("fscoreDfs<-", function(x, ..., value) standardGeneric("fscoreDfs<-"))


#' @export
setGeneric('fscoreDfNames<-', function(x, value)standardGeneric('fscoreDfNames<-'))

#' @export
setMethod('fscoreDfNames', 'SingleCellExperiment', function(x){
    colnames(int_elementMetadata(x)[[.fscore_key]])
}) 


#' @export
setReplaceMethod('fscoreDfNames', c('SingleCellExperiment', 'character'), function(x, value){
    .set_internal_names(x, value,
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.fscore_key,
        unname.key = .unnamed.fscore
    )
})

#' @export
setMethod("fscoreDf", c("SingleCellExperiment", "missing"), function(x, type) {
    .get_internal_missing(x,
        basefun=fscoreDf,
        namefun=fscoreDfNames,
        funstr="fscoreDf"
    )
})

#' @export
setMethod("fscoreDf", c("SingleCellExperiment", "numeric"), function(x, type) {
    out <- .get_internal_numeric(x, type,
        getfun=int_elementMetadata,
        key=.fscore_key,
        funstr="fscoreDf",
        substr="type")

    out
})

#' @export
setMethod("fscoreDf", c("SingleCellExperiment", "character"), function(x, type) {
    out <- .get_internal_character(x, type,
        getfun=int_elementMetadata,
        key=.fscore_key,
        funstr="fscoreDf",
        substr="type",
        namestr="fscoreDfNames")

    out
})


#' @export
setReplaceMethod("fscoreDf", c("SingleCellExperiment", "missing"), function(x, type, ..., value) {
    .set_internal_missing(x, value,
        basefun=`fscoreDf<-`,
        namefun=fscoreDfNames
    )
})


#' @export
setReplaceMethod("fscoreDf", c("SingleCellExperiment", "numeric"), function(x, type, ..., value) {
    .set_internal_numeric(x, type, value,
        getfun = int_elementMetadata,
        setfun = `int_elementMetadata<-`,
        key = .fscore_key,
	convertfun = NULL,
	xdimfun = nrow,
	vdimfun = nrow,
        funstr = "fscoreDf",
	xdimstr = 'nrow',
	vdimstr = 'nrow',
        substr = "type")
})


#' @export
setReplaceMethod("fscoreDf", c("SingleCellExperiment", "character"), function(x, type, ..., value) {
    .set_internal_character(x, type, value,
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.fscore_key,
        convertfun=NULL,
        xdimfun=nrow,
        vdimfun=nrow,
        funstr="fscoreDf",
        xdimstr="nrow",
        vdimstr="nrow",
        substr="type")
})

#' @export
setMethod("fscoreDfs", "SingleCellExperiment", function(x) {
    value <- .get_internal_all(x, 
        getfun=int_elementMetadata, 
        key=.fscore_key)
    value
})


#' @export
setReplaceMethod("fscoreDfs", "SingleCellExperiment", function(x, value) {
    .set_internal_all(x, value,
        getfun=int_elementMetadata,
        setfun=`int_elementMetadata<-`,
        key=.fscore_key,
        convertfun= NULL,
        xdimfun=nrow,
        vdimfun=nrow,
        funstr="fscoreDfs",
        xdimstr="nrow",
        vdimstr="nrow")
})
