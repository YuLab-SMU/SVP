.gsva_key <- 'gsvaExps'

.fscore_key <- 'fscoreDfs'

.sv_key <- 'svDfs'

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
            "({.cls {class(x)}}, type=\"character\", ...)':\n  ",
            "'", index, "' not in '", namestr, "{.cls {class(x)}}.")))
    })
}

.unnamed.gsva <- 'unnamed.gsva'
.unnamed.fscore <- 'unnamed.fscore'
.unnamed.sv <- '.unnamed.sv'

unamekeys <- c(.unnamed.gsva, .unnamed.fscore, .unnamed.sv)
names(unamekeys) <- c(.gsva_key, .fscore_key, .sv_key)

.set_internal_names <- function (x, value, getfun, setfun, key, unname.key=.unnamed.gsva){
    x <- updateObject(x)
    tmp <- getfun(x)
    value <- .clean_internal_names(value, N = ncol(tmp[[key]]),
        msg = "value", unname.key)
    colnames(tmp[[key]]) <- value
    setfun(x, tmp)
}

.unnamed.gsva <- 'unnamed.gsva'

.clean_internal_names <- function(names, N, msg, unname.key){
    if (is.null(names) && N > 0){
        cli::cli_warn(paste0("'", msg, "' is NULL, replacing with '", .unnamed.gsva,"'."))
        names <- paste0(unname.key, seq_len(N))
    }else if (any(empty <- nchar(names) == 0)){
        cli::cli_warn(paste0("'", msg, "' contains empty strings, replacing with '", .unnamed.gsva,"'."))
        names[empty] <- paste0(unname.key, seq_along(sum(empty)))
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

#' @importFrom S4Vectors metadata<- metadata mcols mcols<-
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
            msg = "names(value)", unname.key = unamekeys[key])
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

.check_gsvaexp_columns <- function(main, alt, withDimnames, withColData, 
                                   withSpatialCoords, withImgData, withReducedDim, 
                                   fun = "gsvaExp", vname = "value"){
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
        flag1 <- .check_element_obj(alt, key = 'spatialCoords', basefun = int_colData, namefun = names)
        flag2 <- .check_element_obj(alt, key = 'imgData', basefun = int_metadata, namefun = names)
        if (!withSpatialCoords && flag1 || (!withImgData && flag2)){
            alt <- as(alt, "SingleCellExperiment")
            int_colData(alt)[["spatialCoords"]] <- NULL
            int_metadata(alt)[['imgData']] <- NULL
        }
        if (!withReducedDim){
            reducedDims(alt) <- NULL
        }
    }
    alt
}    

#' @importFrom SpatialExperiment spatialCoords imgData imgData<- spatialCoords<-
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDims reducedDims<-
.fill_gsvaexps_info <- function(out, x, withDimnames, withColData, withSpatialCoords, withImgData, withReducedDim){
    if (withDimnames) {
        colnames(out) <- colnames(x)
    }
    flag1 <- .check_element_obj(x, key = 'spatialCoords', basefun = int_colData, namefun = names)
    if (withColData){
        cmain <- colData(x)
        cout <- colData(out)
        cmain <- cmain[, setdiff(colnames(cmain), colnames(cout)),drop=FALSE]
        if (flag1){
            cmain$sample_id <- NULL
        }
        prep <- cbind(cmain, cout)
        rownames(prep) <- colnames(out)
        colData(out) <- prep
    }
    if (flag1 && withSpatialCoords){
        out$sample_id <- x$sample_id
        out <- as(out, 'SpatialExperiment')
        spatialCoords(out) <- .extract_element_object(x, key = 'spatialCoords', basefun = int_colData, namefun = names)
    }
    flag2 <- .check_element_obj(x, key = 'imgData', basefun = int_metadata, namefun = names)
    if (flag2 && withImgData){
        out <- as(out, 'SpatialExperiment')
        imgData(out) <- .extract_element_object(x, key = 'imgData', basefun = int_metadata, namefun = names)
    }
    if (withReducedDim){
        reducedDims(out) <- reducedDims(x)
    }
    return(out)
}

.sce_to_svpe <- function(sce, gsvaExps = list()){
    if (inherits(sce, "SVPExperiment")){
        return(sce)
    }
    svpe <- new('SVPExperiment', sce)
    int_colData(svpe)[[.gsva_key]] <- new('DFrame', nrows=ncol(svpe))
    gsvaExps(svpe) <- gsvaExps
    return(svpe)
}


.check_element_obj <- function(x, key, basefun, namefun){
    tmp <- basefun(x)
    if (key %in% namefun(tmp)){
        tmp <- tmp[[key]]
        return((inherits(tmp,'matrix') || inherits(tmp, 'DFrame')) && nrow(tmp) > 0)
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


.check_sample_id <- function(x, sampleid){
    sampleid <- unique(sampleid)
    allsample <- .extract_sampleid(x)
    if (is.null(allsample) || sampleid == ".ALLCELL"){
        return(".ALLCELL")
    }
    if (length(sampleid)==1 && sampleid == 'all'){
        return(allsample)
    }
    
    ids <- intersect(sampleid, allsample)

    if (length(ids) < 1){
        cli::cli_abort("The `sample_id` is/are not present in the object.
                        Please check the `sample_id`.", call=NULL)
    }

    if (length(ids) != length(sampleid)){
        cli::cli_inform("Some sample_id are not present in the object.
                        Only using the sample_id, which is/are in the object.")
    }
    return(ids)
}

.extract_sampleid <- function(x){
    unique(colData(x)$sample_id)
}

.tidy_sv_result <- function(x){
    rlang::check_installed(c("tidyr", "tibble"))
    if (length(x) == 1){
        return(x[[1]])
    }

    x <- lapply(x,function(i) i |> tibble::rownames_to_column(var="features")) |>
    dplyr::bind_rows(.id='sample_id') |> 
    dplyr::group_by(.data$features) |> 
    tidyr::nest() |> 
    dplyr::ungroup()
    return(x)
}

.check_features <- function(x, y, prefix){
  x <- unique(x)
  f1 <- match(x, y)
  f1 <- f1[!is.na(f1)]
  if (length(f1) < 1){
      cli::cli_abort(paste0("The `", prefix[1],"` is/are not present in the row names."))
  }
  return(f1)
}


.check_coords <- function(data, 
                          reduction.used, 
                          weight = NULL, 
                          prefix="Or the `weight` should be provided."
  ){
  flag1 <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)

  flag2 <- any(reduction.used %in% reducedDimNames(data))
  coords <- NULL
  if((flag1 || flag2) && is.null(weight)){
      if (flag2){
          coords <- reducedDim(data, reduction.used)
          coords <- coords[,c(1, 2)]
      }
      if (flag1 ){
          coords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
      }
  }else if (all(!flag1, !flag2) && is.null(weight)){
      cli::cli_abort(c("The {.cls {class(data)}} should have 'spatialCoords' or the",
                     paste("reduction result of 'UMAP' or 'TSNE'.", prefix)))
  }

  return(coords)
}

.check_features_in_sce <- function(sce, svg){
  if (inherits(svg, "tbl_df")){
      sce <- sce[svg$features,]
  }else{
      sce <- sce[rownames(svg),]  
  }
  return(sce)
}

.extract.gset <- function(x, gene.name=FALSE){
  if (inherits(x, "list") && !is.null(names(x))){
      return(x)
  }
  if (inherits(x, "GSON")){
      x <- .extract.gset.from.gson(x, gene.name=gene.name)
  }else if (inherits(x, 'character') && file.exists(x)){
      x <- .read.gmt(x)
  }else{
      cli::cli_abort(c("The `gset.idx.list` must be a list which have name (gene set name, such as ",
                       " GO Term name or Reactome Pathway name) or GSON object defined in `gson` package.",
                       "Or the gmt file."))
  }
  return(x)
}


.read.gmt <- function(gmtfile){
    x <- readLines(gmtfile)
    res <- strsplit(x, "\t")
    names(res) <- vapply(res, function(y) y[1], character(1))
    res <- lapply(res, "[", -seq(2))
    return(res)
}

.extract.gset.from.gson <- function(x, gene.name=FALSE){
  gsid2gene <- x@gsid2gene
  gene2name <- x@gene2name
  gsid2name <- x@gsid2name
  nm <- "gsid"
  if (!is.null(gsid2name)){
      y <- dplyr::left_join(gsid2gene, gsid2name, by='gsid')
      nm <- "name"
  }
  gnm <- colnames(gsid2gene)[2]
  if (!is.null(gene2name) && gene.name){
      ind <- colnames(gene2name)[1] |> setNames(gnm)
      y <- dplyr::left_join(y, gene2name, by = c(ind))
      gnm <- colnames(gene2name)[2]
  }else if(gene.name && is.null(gene2name)){
      cli::cli_warn("The `gson` object does not have `SYMBOL` ID.")
  }
  
  x <- y |> 
       dplyr::group_by(!!rlang::sym(nm)) |> 
       dplyr::summarize(.TARGET=list(!!rlang::sym(gnm))) |> 
       dplyr::pull(".TARGET", name=nm)
  return(x)
}
