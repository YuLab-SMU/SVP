#' @export
setGeneric('gsvaExps', function(x,...)standardGeneric('gsvaExps'))

#' @export
setGeneric('gsvaExp', function(x, e, ...)standardGeneric('gsvaExp'))

#' @export
setGeneric('gsvaExpNames', function(x, ...)standardGeneric('gsvaExpNames'))

#' @export
setGeneric('mainGsvaExpName', function(x)standardGeneric('mainGsvaExpName'))

#' @export
setGeneric("gsvaExp<-", function(x, e, withDimnames=TRUE, withColData=FALSE, ..., value) standardGeneric("gsvaExp<-"))

#' @export
setGeneric("gsvaExps<-", function(x, withDimnames=TRUE, withColData=FALSE, ..., value) standardGeneric("gsvaExps<-"))

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

#' @importFrom SingleCellExperiment int_metadata
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
setReplaceMethod('gsvaExp', c('SVPExperiment', 'missing'), function(x, e, withDimnames = TRUE, withColData = FALSE, ..., value){
    .set_internal_missing(x, value,
        withDimnames=withDimnames,
        withColData=withColData,
        basefun=`gsvaExp<-`,
        namefun=gsvaExpNames
    )

})

#' @export
setReplaceMethod('gsvaExp', c('SVPExperiment', 'character'), function(x, e, withDimnames = TRUE, withColData = FALSE, ..., value){
    value <- .check_gsvaexp_columns(x, value, withDimnames=withDimnames, withColData=withColData)
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



#' @export
setReplaceMethod('gsvaExp', c('SVPExperiment', 'numeric'), function(x, e, withDimnames = TRUE, withColData = FALSE, ..., value){
    value <- .check_gsvaexp_columns(x, value, withDimnames=withDimnames, withColData=withColData)

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
setReplaceMethod('gsvaExps', 'SVPExperiment', function(x, withDimnames=TRUE, withColData=FALSE, ..., value){
    if (withDimnames || withColData) {
        for (v in seq_along(value)) {
            value[[v]] <- .check_gsvaexp_columns(x, value[[v]],
                withDimnames=withDimnames, withColData=withColData,
                fun='altExps', vname=sprintf("value[[%s]]", v))
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
setMethod('gsvaExps', 'SVPExperiment', function(x, withDimnames=TRUE, withColData=TRUE, withSpatialCoords=TRUE, withImgData=TRUE, ...){
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
                       withImgData = withImgData)
    }
    return(y)
})

#' @export
setMethod('gsvaExp', c('SVPExperiment', "missing"), function(x, e, withDimnames = TRUE, withColData = TRUE, withSpatialCoords = TRUE, withImgData = TRUE, ...){
    y <- .get_internal_missing(x,
            basefun = gsvaExp,
            namefun = gsvaExpNames,
            funstr = 'gsvaExp',
            withDimnames = withDimnames,
            withColData = withColData,
            withSpatialCoords = withSpatialCoords,
            withImgData = withImgData,
            ...
         )
    return(y)
})

#' @export
setMethod('gsvaExp', c('SVPExperiment', 'numeric'), function(x, e, withDimnames = TRUE, withColData = TRUE, withSpatialCoords = TRUE, withImgData = TRUE, ...){
    y <- .get_internal_numeric(
           x,
           index = e,
           getfun = int_colData,
           key = .gsva_key,
           funstr = 'gsvaExp',
           substr = 'e'
         )
    y <- .get_sce(y)

    y <- .fill_gsvaexps_info(y, x, withDimnames, withColData, withSpatialCoords, withImgData)
    return(y)

})


#' @export
setMethod('gsvaExp', c('SVPExperiment', 'character'), function(x, e, withDimnames = TRUE, withColData = TRUE, withSpatialCoords = TRUE, withImgData = TRUE, ...){
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
    y <- .fill_gsvaexps_info(y, x, withDimnames, withColData, withSpatialCoords, withImgData)
    return(y)
})

