#' @export
setGeneric('plot_point_features', function(x, assayName, features, mapping=NULL, ncol=3, ...)standardGeneric('plot_point_features'))


#' @importFrom ggplot2 aes_string ggplot
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom dplyr left_join
#' @export
setMethod('plot_point_features', c('SpatialExperiment'), function(x, assayName, features, mapping = NULL, ncol = 3){
    x1 <- x[features,]
    x1 <- assay(x1, assayName)
    x1 <- x1 |> as.matrix() |> 
        tibble::as_tibble(rownames='.features') |>
        tidyr::pivot_longer(cols=!".features", 
                            names_to = '.samples', 
                            values_to = '.values'
        )
    coords <- spatialCoords(x) |>
              as.matrix() |>
              tibble::as_tibble(rownames='.samples')

    default_mapping <- aes_string(x=colnames(coords)[3], 
                                  y = colnames(coords)[2], 
                                  color='.values')
    
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
        xlab(NULL)

})

