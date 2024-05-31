#' @title LISAResult
#' @description 
#' Extracting the result of `runLISA()`
#' @param x object \linkS4class{SingleCellExperiment}.
#' @param type character, the name of \code{method} parameter of \code{runLISA()} combining 
#' \code{.SVP},so it can be one of \code{localG.SVP} or \code{localmoran.SVP},
#' default is NULL.
#' @param features character or index which have been specified in \code{features} of 
#' `code{runLISA()}` and \code{action='add'}, default is NULL.
#' @param ... additional parameter, meaningless now.
#' @return a data.frame or SimpleList.
#' @export
LISAResult <- function(x, type = NULL, features=NULL, ...){
    if (!inherits(x, "SingleCellExperiment")){
        cli::cli_abort("The `x` should be a `SingleCellExperiment` object!")
    }
    if (!'localResults' %in% names(x@int_colData)){
        cli::cli_abort(c("The `x` (SingleCellExperiment) should have preformed `runLISA()` function", 
                         "with `action = 'add'`."))
    }
    x <- x@int_colData$localResults |> 
          lapply(function(x) x|> as.list() |> SimpleList())  

    if (length(type) > 2){
        cli::cli_warn("The `type` should be a one length vector. The first element is used!")
        type <- type[1]
    }
    
    if (is.null(type)){
        type <- names(x)[1]
        cli::cli_warn(paste0("The `type` are not specified. The first `type='", 
                              type,"'` are extracted!"))
    }
    
    x <- x[[type]]

    if (is.null(features)){
        cli::cli_warn("The `features` are not specified. All features are extracted!")
        features <- names(x)
    }
    if (is.numeric(features)){
        features <- names(x)[features]
    }
    features <- intersect(features, names(x))
    if (length(features) > 1){
        return(x[features])
    }
    x[[features]]
}

#' convert the square matrix to long tidy table
#' @description
#' This function is designed to convert the output of \code{runGLOBALBV},
#' \code{fast_cor} or the matrix output of \code{cor} to long tidy table.
#' @param x list or matrix object, which is the output of \code{runGLOBALBV},
#' \code{fast_cor} or the matrix output of \code{cor}.
#' @param listn list object, which must have name, and the element must
#' from the row names of \code{x} or \code{x[[1]]} (when \code{x} is a list)
#' default is NULL.
#' @return a long tidy table
#' @export
#' @examples
#' example(fast_cor, echo=FALSE)
#' x <- as_tbl_df(res)
#' head(x)
#' x2 <- as_tbl_df(res2)
#' head(x2)
as_tbl_df <- function(x, listn=NULL){
    if (inherits(x, 'list') && length(x)==2){
        rmat <- x[[1]]
        pval <- x[[2]]
        nm <- names(x) 
    }else{
        rmat <- x
        pval <- NULL
        nm <- c('val')
    }
    rmat <- .internal.as_tbl_df(rmat)
    colnames(rmat) <- c("x", "y", nm[1])
    if (!is.null(pval)){
        pval <- .internal.as_tbl_df(pval)
        colnames(pval) <- c("x", "y", nm[2])
        rmat <- dplyr::left_join(rmat, pval, by=c("x", "y"))
    }
    if (!is.null(listn)){
        if (inherits(listn, "list") && !is.null(names(listn))){
            rmat <- .add_list(rmat, listn)
        }else{
            cli::cli_warn("The `listn` must be a list object with name.")
        }
    }
    rmat <- tibble::as_tibble(rmat)
    return(rmat)
}

