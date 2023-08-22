#' @export
setGeneric('plot_point_features', function(x, 
                                           assayName, 
                                           features, 
                                           reducedDim.type = NULL,
                                           dims = c(1, 2),
                                           mapping=NULL, 
                                           ncol=3, 
                                           ...)
  standardGeneric('plot_point_features')
)


#' @importFrom ggplot2 aes_string ggplot
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom dplyr left_join
#' @importFrom SummarizedExperiment rowData assayNames colData
#' @importFrom SingleCellExperiment reducedDim
#' @export
setMethod('plot_point_features', c('SingleCellExperiment'), 
          function(
            x, 
            assayName, 
            features, 
            reducedDim.type=NULL,
            dims = c(1, 2), 
            mapping = NULL, 
            ncol = 3, 
            ...){
    x1 <- x[features,]
    x1 <- assay(x1, assayName)
    if (is.numeric(assayName) || is.logical(assayName)){
        assayName <- assayNames(x)[assayName]
    }
    x1 <- x1 |> as.matrix() |> 
        tibble::as_tibble(rownames='.features') |>
        tidyr::pivot_longer(cols=!".features", 
                            names_to = '.samples', 
                            values_to = assayName
        )
    flag1 <- .check_element_obj(x, key='spatialCoords', basefun = int_colData, namefun = names) 
    if (flag1 && is.null(reducedDim.type)){
        coords <- .extract_element_object(x, key = 'spatialCoords', basefun = int_colData, namefun = names) 
        coords.obj <- coord_fixed()
    }else{
        coords <- reducedDim(x, reducedDim.type)
        colnames(coords) <- paste0(reducedDim.type, seq_len(ncol(coords)))
        if (!is.null(dims)){
            if (length(dims)>1){
                coords <- coords[, dims[seq_len(2)]]
            }else{
                cli::cli_abort(c('The {.arg dims} should be a two-length or larger two-length vector.'))
            }
        }
        coords.obj <- NULL
    }
    coords <- as.matrix(coords) |> 
              tibble::as_tibble(rownames = '.samples')

    default_mapping <- aes_string(x=colnames(coords)[2], 
                                  y = colnames(coords)[3], 
                                  color= assayName)
    
    feature.da <- rowData(x) |> 
                  as.matrix() |> 
                  tibble::as_tibble(rownames='.features')
    sample.da <- colData(x) |>
                 as.matrix() |>
                 tibble::as_tibble(rownames = '.samples')
    x1 <- x1 |>
          dplyr::left_join(coords, by = '.samples') |> 
          dplyr::left_join(feature.da, by = '.features') |>
          dplyr::left_join(sample.da, by = '.samples')

    if (!is.null(mapping)){
        mapping <- modifyList(default_mapping, mapping)
    }else{
        mapping <- default_mapping
    }

    ggplot(x1, mapping = mapping) +
        geom_point() +
        .feature_setting(features, ncol) +
        ylab(NULL) +
        xlab(NULL) +
        coords.obj +
        scale_color_gradient2()

})


#' @export
setMethod('plot_point_features', 'SpatialExperiment', 
          function(
            x,
            assayName,
            features,
            reducedDim.type=NULL,
            dims = c(1, 2),
            mapping = NULL,
            ncol = 3,
            ...){
    callNextMethod()
})
