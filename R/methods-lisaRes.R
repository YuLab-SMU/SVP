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
#' @importFrom S4Vectors SimpleList
#' @export
#' @examples
#' data(hpda_spe_cell_dec)
#' hpda_spe_cell_dec <- hpda_spe_cell_dec |> 
#'                      runLISA(features = 'Cancer clone A', 
#'                              assay.type = 'affi.score', 
#'                              method = 'localG',
#'                              action = 'add'
#'                      ) 
#' hpda_spe_cell_dec <- hpda_spe_cell_dec |>
#'                      runLISA(features = 'Cancer clone A', 
#'                              assay.type = 'affi.score', 
#'                              method = 'localmoran',
#'                              action = 'add'
#'                      )
#' local.G <- LISAResult(hpda_spe_cell_dec, 
#'                       type='localG.SVP', features='Cancer clone A'
#'            )
#' localmoran <- LISAResult(hpda_spe_cell_dec, 
#'                          type = 'logcalmoran.SVP', 
#'                          features = 'Cancer clone A'
#'            )
#' hpda_spe_cell_dec |> LISAResult() |> head()
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


#' @title convert LISA result to SVPExperiment.
#' @description
#' convert the Gi for \code{runLISA} result or LocalLee for \code{runLOCALBV} result 
#' to a \linkS4class{SVPExperiment}.
#' @rdname LISAsce-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object.
#' @param lisa.res list the result of \code{runLISA} or \code{runLOCALBV}.
#' @param gsvaexp.name character the name of gsveExp for the LISA result, default is "LISA".
#' @param ... currently meaningless.
#' @seealso [`runLISA`] and [`runLOCALBV`]
#' @examples
#' data(hpda_spe_cell_dec)
#' lisa.res12 <- hpda_spe_cell_dec |>
#'    runLISA(
#'      features = c(1, 2, 3),
#'      assay.type = 'affi.score',
#'      weight.method = "knn",
#'      k = 10,
#'      action = 'get',
#'    )
#' hpda_spe_cell_dec <- LISAsce(hpda_spe_cell_dec, lisa.res12)
#' hpda_spe_cell_dec
#' gsvaExp(hpda_spe_cell_dec, 'LISA')
#' localbv.res1 <- hpda_spe_cell_dec |> runLOCALBV(
#'           features1 = 'Cancer clone A',
#'           features2 = 'Cancer clone B',
#'           assay.type='affi.score'
#'         )
#' hpda_spe_cell_dec <- LISAsce(hpda_spe_cell_dec, localbv.res1, 'LOCALBV')
#' gsvaExp(hpda_spe_cell_dec, 'LOCALBV')
setGeneric('LISAsce',
  function(
    data,
    lisa.res,
    gsvaexp.name = "LISA",
    ...
  )
  standardGeneric('LISAsce')
)

#' @rdname LISAsce-method
#' @aliases LISAsce,SingleCellExperiment
#' @export LISAsce
setMethod("LISAsce", "SingleCellExperiment", 
  function(
    data, 
    lisa.res, 
    gsvaexp.name = "LISA", 
    ...
  ){
    if (missing(lisa.res)){
      if (is.null(data@int_colData$localResults)){
         cli::cli_warn(c("The {.var lisa.res} should be provided, when the {.var data} do not contain ",
                         "the result of runLISA or runLOCALBV (speficied action='add')."))
      }
      lisa.res <- data@int_colData$localResults |> lapply(function(x) x|> as.list() |> SimpleList())
      lisa.res <- lisa.res[[1]]
    }
    flag <- .check_localbv_res(lisa.res)
    dat <- .build_sce_assays(lisa.res, flag)
    data <- .sce_to_svpe(data)
    gsvaExp(data, gsvaexp.name) <- SingleCellExperiment(assays = dat)
    return(data)
})

#' convert the square matrix to long tidy table
#' @description
#' This function is designed to convert the output of \code{runGLOBALBV},
#' \code{fast_cor} or the matrix output of \code{cor} to long tidy table.
#' @param x list or matrix object, which is the output of \code{runGLOBALBV},
#' \code{fast_cor} or the matrix output of \code{cor}.
#' @param listn list object, which must have name, and the element must
#' from the row names of \code{x} or \code{x[[1]]} (when \code{x} is a list)
#' default is NULL.
#' @param diag logical whether include the diagonal (only work when the cor 
#' matrix is square), default is TRUE.
#' @param rmrd logical whether remove of redundancy when the correlation matrix
#' is a square matrix, default is TRUE.
#' @param flag.clust logical whether perform the hierarchical cluster analysis
#' to obtain the label for visualization.
#' @param dist.method the distance measure to be used, only work when \code{flag.clust = TRUE}.
#' It must be one of \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"},
#' \code{"binary"} or \code{"minkowski"}.
#' @param hclust.method the agglomeration method to be used, only work with \code{flag.clust=TRUE}.
#' This should be (an unambiguous abbreviation of) one of \code{"ward.D"}, \code{"ward.D2"}, 
#' \code{"single"}, \code{"complete"}, \code{"average"} (= UPGMA), \code{"mcquitty"} (= WPGMA), 
#' \code{"median"} (= WPGMC) or \code{"centroid"} (= UPGMC).
#' @return a long tidy table
#' @export
#' @importFrom stats dist hclust
#' @examples
#' library(ggplot2)
#' library(ggtree)
#' library(aplot)
#' example(fast_cor, echo=FALSE)
#' x <- as_tbl_df(res)
#' head(x)
#' xx <- as_tbl_df(res, flag.clust = TRUE, 
#'                 dist.method = 'euclidean', hclust.method = 'average')
#' p1 <- ggplot(xx, mapping = aes(x=x,y=y,color=r,size=abs(r))) +
#'       geom_point() + xlab(NULL) + ylab(NULL) +
#'       guides(y=guide_axis(position='right'))
#' p2 <- res$r |> dist() |> hclust(method = 'average') |> 
#'       ggtree(layout='den', branch.length='none', ladderize=FALSE)
#' p3 <- res$r |> t() |> dist() |> hclust(method = 'average') |>
#'       ggtree(branch.length = 'none', ladderize = FALSE)
#' p4 <- p1 |> insert_left(p3, width=.12) |> insert_top(p2, height=.12)
#' aplot::plot_list(p1, p4)
#' x2 <- as_tbl_df(res2)
#' head(x2)
#' f1 <- ggplot(x2, aes(x=x, y=y, color=r, size=abs(r))) + geom_point() +
#'       xlab(NULL) + ylab(NULL) + 
#'       guides(x=guide_axis(position='top', angle=45), 
#'              y=guide_axis(position='right'))
#' f2 <- res2$r |> t() |> dist() |> hclust(method = 'average') |>
#'       ggtree(branch.length = 'none', ladderize=FALSE)
#' f3 <- f1 |> aplot::insert_left(f2, width=.12)
#' xx2 <- as_tbl_df(res2, 
#'                 flag.clust = TRUE,
#'                 dist.method = 'euclidean', 
#'                 hclust.method = 'average'
#'       )
#' ff1 <- ggplot(xx2, mapping = aes(x=x,y=y, color=r,size=abs(r))) +
#'        geom_point() + xlab(NULL) + ylab(NULL) +
#'        guides(x=guide_axis(position='top', angle=45),
#'               y=guide_axis(position='right'))
#' ff3 <- ff1 |> aplot::insert_left(f2, width = .12)
#' aplot::plot_list(f3, ff3)
as_tbl_df <- function(x, 
                      listn = NULL, 
                      diag = TRUE, 
                      rmrd = TRUE,
                      flag.clust = FALSE, 
                      dist.method = 'euclidean',
                      hclust.method = 'average'
                     ){
    if (inherits(x, 'list') && length(x)==2){
        rmat <- x[[1]]
        pval <- x[[2]]
        nm <- names(x) 
    }else{
        rmat <- x
        pval <- NULL
        nm <- c('val')
    }
    rmat <- .internal.as_tbl_df(rmat, diag = diag, rmrd = rmrd, flag.clust=flag.clust, 
                                dist.method=dist.method, hclust.method=hclust.method)
    colnames(rmat) <- c("x", "y", nm[1])
    if (!is.null(pval)){
        pval <- pval[levels(rmat$x), levels(rmat$y), drop=FALSE]
        pval <- .internal.as_tbl_df(pval, diag = diag, rmrd = rmrd)
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

.build_sce_assays <- function(x, flag){
    ind <- seq(1)
    nm <- "Gi"
    if (flag){
       ind <- seq(2)
       nm <- c("LocalLee", "Gi")
    }
    x <- lapply(ind, function(i).build_matrix(x, ind = i))
    names(x) <- nm
    return(x)
}

.build_matrix <- function(x, ind=1){
    nm <- names(x)
    x <- lapply(x, function(i)i[,ind,drop=FALSE])
    x <- do.call('cbind', x) |> setNames(nm) |> t()
    x <- Matrix::Matrix(x, sparse=TRUE)
    return(x)   
}

.check_localbv_res <- function(x){
    "LocalLee" %in% colnames(x[[1]])
}
