#' @title runCORR
#' @description
#' This function to perform the correlation of the features in main experiment or 
#' features of gsva experiment.
#' @rdname runCORR-method
#' @param data a \linkS4class{SingleCellExperiment} object with contains \code{UMAP} or \code{TSNE},
#' or a \linkS4class{SpatialExperiment} object, or a \linkS4class{SVPExperiment} object with specified
#' \code{gsvaexp} argument.
#' @param features1 the features name data object (only supporting character), default is NULL, 
#' see also \code{features2} parameter.
#' @param features2 character, if \code{features1} is not NULL, and \code{features2} is NULL, 
#' only the \code{features1} are analyzed, if \code{features1} is NULL, and \code{features2} is
#' is not NULL, the \code{features2} are analyzed, if \code{features2} is also NULL, all of features in the
#' \code{data} object will be analyzed. If \code{features2} and \code{features1} are not NULL, the bivariate
#' spatial autocorrelation analysis will be performed between the \code{features1} and \code{features2}.
#' default is \code{NULL}.
#' @param assay.type which expressed data to be pulled to run, default is \code{logcounts}.
#' @param method character should be one of the \code{spearman}, \code{pearson} and \code{bicorr}, default is \code{'spearman'}.
#' @param alternative indicates the alternative hypothesis and must be one of \code{"two.sided"}, \code{"greater"} or 
#' \code{"less"}. You can specify just the initial letter. \code{"greater"} corresponds to positive association, 
#' \code{"less"} to negative association, default is \code{"two.sided"}.
#' @param add.pvalue logical whether calculate the pvalue, which is calculated with permutation test. So it might
#' be slow, default is \code{FALSE}, which the pvalue of result will be NULL.
#' @param action character, which should be one of \code{'only'} and \code{'get'}, default is \code{"get"}.
#' If \code{action='only'}, it will return a long tidy table contains the correlation for each feature pairs.
#' If \code{action='get'}, it will return a list containing the correlation matrix and pvalue matrix (if \code{add.pvalue=TRUE}).
#' @param verbose logical whether print the help information, default is TRUE.
#' @param gsvaexp character the one character from the name of `gsvaExpNames(data)`, default is NULL. If \code{data}
#' is \linkS4class{SVPExperiment}, and the parameter is specified simultaneously. the `features` (Usually genes)
#' from the displayed class, and `gsvaexp.features` from name in `rownames(gsvaExp(data, gsvaexp))` will be
#' performed the analysis.
#' @param gsvaexp.assay.type character the assay name in the `assays(gsvaExp(data, gsvaexp))`, default is NULL, 
#' which works with \code{gsvaexp} parameter.
#' @param gsvaexp.features character the name from the `rownames(gsvaExp(data, gsvaexp))`. If \code{gsvaexp} is
#' specified and \code{data} is \linkS4class{SVPExperiment}, it should be provided. Default is NULL.
#' @param across.gsvaexp logical whether only calculate the relationship of features between the multiple `gsvaExps`
#' not the internal features of gsvaExp. For example, \code{'a'} and \code{'b'} features are from the \code{'AB'} `gsvaExp`,
#' \code{'c'} and \code{'d'} features are from the \code{'CD'} `gsvaExp`. When \code{across.gsvaexp=TRUE} and 
#' \code{gsvaexp.features = c('a', 'b', 'c', 'd')} and \code{gsvaexp = c('AB', 'CD')}, Only the relationship of
#' \code{a} and \code{c}, \code{a} and \code{d}, \code{b} and \code{c}, and \code{b} and \code{d} will be calculated.
#' default is TRUE.
#' @param ... additional parameters the parameters which are from the weight.method function.
#' @return long tidy table or list see also the help information of \code{action} argument.
#' @seealso [`runCORR`] to explore the global bivariate relationship in the spatial space.
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(hpda_spe_cell_dec)
#' rownames(hpda_spe_cell_dec) |> head()
#' res1 <- runCORR(hpda_spe_cell_dec, 
#'                 features1 = "Ductal APOL1 high-hypoxic", 
#'                 features2 = c('Cancer clone A', "Cancer clone B"), 
#'                 assay.type = 'affi.score', 
#'                 action='only'
#'         )
#' res1
#' res2 <- runCORR(hpda_spe_cell_dec, 
#'                 features1 = c("Acinar cells", 
#'                               "Ductal APOL1 high-hypoxic", 
#'                               "Cancer clone A",
#'                               "Cancer clone B"), 
#'                 assay.type = 1, 
#'                 action = 'get'
#'         )
#' res2
setGeneric('runCORR',
  function(
    data,
    features1 = NULL,
    features2 = NULL,
    assay.type = 'logcounts',
    method = c("spearman", "pearson", "bicorr"),
    alternative = c('greater', 'two.sided', 'less'),
    add.pvalue = FALSE,
    action = c("get", "only"),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    across.gsvaexp = TRUE,
    ...
  )
  standardGeneric('runCORR')
)

#' @rdname runCORR-method
#' @aliases runCORR,SingleCellExperiment
#' @export runCORR
setMethod("runCORR", "SingleCellExperiment", function(
    data, 
    features1 = NULL,
    features2 = NULL,
    assay.type = "logcounts",
    method = c("spearman", "pearson", "bicorr"),
    alternative = c("greater", "two.sided", "less"),
    add.pvalue = FALSE,
    action = c("get", "only"),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    across.gsvaexp = TRUE,
    ...
  ){
  
  method <- match.arg(method)
  action <- match.arg(action)
  alternative <- match.arg(alternative)

  if (is.null(assay.type)){
      assay.type <- 1
  }

  x <- assay(data, assay.type)
  features <- union(features1, features2)

  if (length(features) >= 1){
      x <- x[features, ,drop=FALSE]
  }
  
  res <- .internal.runCORR(x, features1, features2, method,
                           alternative, add.pvalue, NULL, across.gsvaexp)
  if (action == 'only'){
      res <- as_tbl_df(res)
  }
  return(res)
})


#' @rdname runCORR-method
#' @aliases runCORR,SVPExperiment
#' @export runCORR
setMethod("runCORR", "SVPExperiment", function(
    data,
    features1 = NULL,
    features2 = NULL,
    assay.type = "logcounts",
    method = c("spearman", "pearson", "bicorr"),
    alternative = c('greater', 'two.sided', 'less'),
    add.pvalue = FALSE,
    action = c("get", "only"),
    verbose = TRUE,
    gsvaexp = NULL,
    gsvaexp.assay.type = NULL,
    gsvaexp.features = NULL,
    across.gsvaexp = TRUE,
    ...
  ){

    if (!is.null(gsvaexp)){
       if (verbose){
          cli::cli_inform(c("The {.var gsvaexp} is specified, the specified {.var gsvaExp} will be",
                           "used to perform the analysis."))
       }
       method <- match.arg(method)
       action <- match.arg(action)
       alternative <- match.arg(alternative)

       if (is.null(assay.type)){
           assay.type <- 1
       }

       x <- assay(data, assay.type)
       if (is.null(features1) && is.null(features2) && is.null(gsvaexp.features)){
           cli::cli_abort(c("The {.var gsvaexp} is specified, and the `data` is {.cls {class(data)}}.",
                            "The `features1`, `features2` and `gsvaexp.features` should not be `NULL` simultaneously."))
       }
       if (!is.null(features1) && !is.null(features2)){
           cli::cli_inform(c("The `data` is {.cls {class(data)}} class, the `features1` and `features2` which are from ",
                            "the displayed object, will be merged to perform analysis with `gsvaexp.features`."))
       }
       if (is.null(gsvaexp.features) && !is.null(gsvaexp)){
           cli::cli_abort(c(paste0("The `data` is {.cls {class(data)}} class and `gsvaexp='", gsvaexp,"'`."),
                          " The `gsvaexp.features` should be provided, you can use ",
                          "`rownames(gsvaExp(data)) |> head()` to view them."))
       }

       features1 <- union(features1, features2)
       features2 <- gsvaexp.features

       if (length(features1) >= 1){
           x <- x[features1, , drop=FALSE]
       }

       x2 <- .extract_gsvaExp_assay(data, gsvaexp, gsvaexp.assay.type)
       x2 <- x2[features2, , drop=FALSE]
       x <- rbind(x, x2)

       listn <- .generate_feature_listn(data, features1, features2, gsvaexp)
       res <- .internal.runCORR(x, features1, features2, method,
                                alternative, add.pvalue, listn, across.gsvaexp)
       if (action == 'only'){
           res <- as_tbl_df(res, listn)
       }
    }else{
       res <- callNextMethod()
    }
    return(res)
})
