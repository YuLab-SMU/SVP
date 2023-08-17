.gsva_key <- 'gsvaExps'

#' @importFrom BiocGenerics updateObject 
# refering to the internal functions of SingleCellExperiment
.get_internal_all <- function(x, getfun, key, to.class='SimpleList'){
    x <- updateObject(x)
    as(getfun(x)[[key]], to.class)
}

.get_sce <- function(x){
    x@sce
}

setMethod("length", "SCEByColumn", function(x) ncol(.get_sce(x)))

#' @importFrom methods initialize
setMethod("[", "SCEByColumn", function(x, i, j, ..., drop=FALSE) {
    initialize(x, sce=.get_sce(x)[,i])
})

setReplaceMethod("[", "SCEByColumn", function(x, i, j, ..., value) {
    left <- .get_sce(x)
    left[,i] <- .get_sce(value)
    initialize(x, sce=left)
})

setMethod("c", "SCEByColumn", function(x, ...) {
    gathered <- lapply(list(x, ...), .get_sce)
    initialize(x, sce=do.call(cbind, gathered))
})

setMethod("names", "SCEByColumn", function(x) colnames(.get_sce(x)))

setReplaceMethod("names", "SCEByColumn", function(x, value) {
    colnames(x@sce) <- value
    x
})

#' @importFrom cli cli_abort
.get_internal_missing <- function (x, basefun, namefun, funstr, ...){
    if (!length(namefun(x)) > 0){
        cli_abort(c(paste0("no available entries for '", funstr, "({.cls {as.character(class(x))}} ...)'")))
    }
    basefun(x, 1L, ...)
}

.get_internal_numeric <- function(x, index, getfun, key, funstr, substr){
    x <- updateObject(x)
    internals <- getfun(x)[[key]]
    tryCatch({
        internals[, index]
    }, error = function(err) {
        cli::cli_abort(c(paste0("invalid subscript '", substr, "' in '", funstr,
             "({.cls {class(x)}}, type=\"numeric\", ...)':\n  ",
            conditionMessage(err))))
    })
}

.get_internal_character <- function (x, index, getfun, key, funstr, substr, namestr){
    x <- updateObject(x)
    internals <- getfun(x)[[key]]
    tryCatch({
        internals[, index]
    }, error = function(err) {
        cli::cli_abort(c(paste0("invalid subscript '", substr, "' in '", funstr,
            "{.cls {class(x)}}, type=\"character\", ...)':\n  ",
            "'", index, "' not in '", namestr, "{.cls {class(x)}}.")))
    })
}


.set_internal_names <- function (x, value, getfun, setfun, key){
    x <- updateObject(x)
    tmp <- getfun(x)
    value <- .clean_internal_names(value, N = ncol(tmp[[key]]),
        msg = "value")
    colnames(tmp[[key]]) <- value
    setfun(x, tmp)
}

.unnamed.gsva <- 'unnamed.gsva'

.clean_internal_names <- function(names, N, msg){
    if (is.null(names) && N > 0){
        cli::cli_warn(paste0("'", msg, "' is NULL, replacing with '", .unnamed.gsva,"'."))
        names <- paste0(.unnamed.gsva, seq_len(N))
    }else if (any(empty <- nchar(names) == 0)){
        cli::cli_warn(paste0("'", msg, "' contains empty strings, replacing with '", .unnamed.gsva,"'."))
        names[empty] <- paste0(.unnamed.gsva, seq_along(sum(empty)))
    }
    names
}

.set_internal_missing <- function(x, value, ..., basefun, namefun){
    if (length(namefun(x))){
        type <- 1L
    }else{
        type <- paste0(.unnamed.gsva, 1L)
    }
    basefun(x, type, ..., value = value)
}

SCEByColumn <- function(sce)new('SCEByColumn', sce = sce)

.set_internal_character <- function (x, type, value, getfun, setfun, key, convertfun, xdimfun,
    vdimfun, funstr, xdimstr, vdimstr, substr){
    x <- updateObject(x)
    if (!is.null(value)){
        if (!is.null(convertfun)){
            value <- convertfun(value)
        }
        if (!identical(vdimfun(value), xdimfun(x))) {
            cli::cli_abort(c(paste0("invalid 'value' in '", funstr, "({.cls {class(x)}}, type=\"character\") <- value'.",
                           "'value' should have number of ", vdimstr, " equal to '", xdimstr, "(x)'")))
        }
    }
    internals <- getfun(x)
    internals[[key]][[type]] <- value
    setfun(x, internals)
}

.set_internal_numeric <- function (x, type, value, getfun, setfun, key, convertfun, xdimfun, vdimfun, funstr, xdimstr, vdimstr, substr){
    x <- updateObject(x)
    if (!is.null(value)){
        if (!is.null(convertfun)) {
            value <- convertfun(value)
        }
        if (!identical(vdimfun(value), xdimfun(x))) {
            cli::cli_abort(paste0("invalid 'value' in '", funstr, "({.cls {class(x)}}, type=\"numeric\") <- value'.", 
                           "'value' should have number of ", vdimstr, " equal to '", xdimstr, "(x)'"))
        }
    }
    internals <- getfun(x)
    if (type[1] > ncol(internals[[key]])){
        cli::cli_abort(c(paste0("'", substr, "' out of bounds in '", funstr,"({.cls {class(x)}}, type='numeric'")))
    }
    internals[[key]][[type]] <- value
    setfun(x, internals)
}    

.set_internal_all <- function(x, value, getfun, setfun, key, 
                              convertfun, xdimfun, vdimfun, 
                              funstr, xdimstr, vdimstr){
    x <- updateObject(x)
    if (length(value) == 0L) {
        collected <- getfun(x)[, 0]
    }
    else {
        original <- value
        if (!is.null(convertfun)) {
            value <- lapply(value, convertfun)
        }
        N <- vapply(value, vdimfun, 0L)
        if (!all(N == xdimfun(x))) {
            cli::cli_abort(c(paste0("invalid 'value' in '", funstr, "(", "{.cls {class(x)}} )<- value'.", 
                             "each element of 'value' should have number of ", vdimstr, " equal to '", xdimstr, "(x)'")))
            
        }
        names(value) <- .clean_internal_names(names(value), N = length(value),
            msg = "names(value)")
        collected <- do.call(DataFrame, c(lapply(value, I), list(row.names = NULL,
            check.names = FALSE)))
        if (is(original, "Annotated")) {
            metadata(collected) <- metadata(original)
        }
        if (is(original, "Vector")) {
            mcols(collected) <- mcols(original)
        }
    }
    tmp <- getfun(x)
    tmp[[key]] <- collected
    setfun(x, tmp)
}

.check_gsvaexp_columns <- function(main, alt, withDimnames, withColData, fun = "gsvaExp", vname = "value"){
    if (!is.null(alt)){
        if (withDimnames) {
            if (!identical(colnames(main), colnames(alt))) {
                msg <- paste0("'colnames(", vname, ")' are not the same as 'colnames(x)' for '",
                  fun, "<-'. This will be an error in the next release of Bioconductor.")
                cli::cli_warn(msg)
            }
        }
        if (withColData) {
            main.cd <- colData(main)
            ncd <- ncol(main.cd)
            alt.cd <- colData(alt)
            acd <- ncol(alt.cd)
            keep <- seq_len(acd) <= ncd
            if (acd < ncd || !identical(alt.cd[, keep, drop = FALSE], main.cd)){
                cli::cli_warn(paste0("left-most columns of 'colData(", vname,
                  ")' should be the same as 'colData(x)' when 'withColData=TRUE'"))
            }else{
                colData(alt) <- alt.cd[, !keep, drop = FALSE]
            }
        }
    }
    alt
}    

#' @importFrom SpatialExperiment spatialCoords imgData
.fill_gsvaexps_info <- function(out, x, withDimnames, withColData, withSpatialCoords, withImgData){
    if (withDimnames) {
        colnames(out) <- colnames(x)
    }
    if (withColData) {
        prep <- cbind(colData(x), colData(out))
        rownames(prep) <- colnames(out)
        colData(out) <- prep
    }
    flag1 <- .check_element_obj(x, key = 'spatialCoords', basefun = int_colData, namefun = names) 
    if (flag1 && withSpatialCoords){
        out <- as(out, 'SpatialExperiment')
        spatialCoords(out) <- .extract_element_object(x, key = 'spatialCoords', basefun = int_colData, namefun = names)
    }
    flag2 <- .check_element_obj(x, key = 'imgData', basefun = int_metadata, namefun = names)
    if (flag2 && withImgData){
        out <- as(out, 'SpatialExperiment')
        imgData(out) <- .extract_element_object(x, key = 'imgData', basefun = int_metadata, namefun = names)
    }
    return(out)
}

.sce_to_svpe <- function(sce, gsvaExps = list()){
    svpe <- new('SVPExperiment', sce)
    int_colData(svpe)[[.gsva_key]] <- new('DFrame', nrows=ncol(svpe))
    gsvaExps(svpe) <- gsvaExps
    return(svpe)
}

#.check_element_obj(x, key='spatialCoords', basefun=int_colData, namefun = names)

.check_element_obj <- function(x, key, basefun, namefun){
    tmp <- .extract_element_object(x=x, basefun = basefun, namefun = namefun)
    if (key %in% namefun(tmp)){
        tmp <- tmp[[key]]
        return(inherits(tmp,'matrix') && !all(is.na(tmp)))
    }else{
        return(FALSE)
    }
}

.extract_element_object <- function(x, key, basefun, namefun){
    tmp <- basefun(x)
    if (!missing(key) && key %in% namefun(tmp)){
        return(tmp[[key]])
    }else{
        return(tmp)
    }
}

#' @importFrom rlang check_installed
.run_sv <- function(x, svgfun, ...){
    params <- list(...)
    params <- c(list(x), params)
    if (missing(svgfun)){
        check_installed('nnSVG', "for svp().", action=BiocManager::install)
        svgfun <- nnSVG::nnSVG
        if ('threads' %in% names(params)){
            names(params)[names(params)=='threads'] <- 'n_threads'
        }
        names(params) <- gsub('sv.', '', names(params))
        res <- do.call(svgfun, params)
        return(res)
    }
}

.tidy_res.sv <- function(x, y){
    x <- S4Vectors::merge(x, y, by = 0, all = TRUE)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
    return(x)
}
