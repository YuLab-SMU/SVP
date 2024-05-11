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
