#' @export
setGeneric('plot_point_features', function(x, 
                                           assayName, 
                                           features = NULL, 
                                           reducedDim.type = NULL,
                                           dims = c(1, 2),
                                           mapping=NULL, 
                                           ncol = 3, 
                                           image.plot = TRUE,
                                           image.rotate.degree = NULL,
                                           image.mirror.axis = 'h',
                                           remove.point = FALSE,
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
            features = NULL, 
            reducedDim.type=NULL,
            dims = c(1, 2), 
            mapping = NULL, 
            ncol = 3,
            image.plot = TRUE,
            image.rotate.degree = NULL,
            image.mirror.axis = 'h',
            remove.point = FALSE,
            ...){

    flag1 <- .check_element_obj(x, key='spatialCoords', basefun = int_colData, namefun = names)
    flag2 <- .check_element_obj(x, key= 'imgData', basefun = int_metadata, namefun = names)

    if (flag1 && is.null(reducedDim.type) && (!image.plot || !flag2)){
        coords <- .extract_element_object(x, key = 'spatialCoords', basefun = int_colData, namefun = names) 
        coords.obj <- coord_fixed()
    }else if (flag2 && image.plot && is.null(reducedDim.type)){
        coords <- .extract_element_object(x, key = 'spatialCoords', basefun = int_colData, namefun = names)
        coords <- .adjust_coords_by_image(coords, x, ...)
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

    feature.da <- rowData(x) |>
                  as.matrix() |>
                  tibble::as_tibble(rownames='.features')
    sample.da <- colData(x) |>
                 as.matrix() |>
                 tibble::as_tibble(rownames = '.samples')

    default_mapping <- aes_string(x = colnames(coords)[2], 
                                  y = colnames(coords)[3])
    
    if (is.null(features) || missing(features)){
        x1 <- coords
    }else{
        x1 <- x[features,]
        if (is.numeric(features)){
            features <- rownames(x1)
        }
        if (missing(assayName)){assayName <- 1L}
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

        x1 <- x1 |>
              dplyr::left_join(coords, by = '.samples') |>
              dplyr::left_join(feature.da,  by = '.features')

        x1$`.features` <- factor(x1$`.features`, levels = features)

        default_mapping <- modifyList(default_mapping, aes_string(color = assayName))
    }

    x1 <- x1 |>
          dplyr::left_join(sample.da, by = '.samples')

    if (!is.null(mapping)){
        mapping <- modifyList(default_mapping, mapping)
    }else{
        mapping <- default_mapping
    }

    p <- ggplot(x1, mapping = mapping) 
    
    if (flag2 && image.plot && is.null(reducedDim.type)){
        annot.img <- .build_image_annotation(x, 
                                             rotate.degree = image.rotate.degree, 
                                             mirror.axis = image.mirror.axis, 
                                             ...)
        p <- p + annot.img
    }

    if (!remove.point || !(is.null(features) || missing(features))){
        p <- p + geom_point()
    }else{
        p <- p + geom_blank()
    }

    p <- p + 
        .feature_setting(features, ncol) +
         ylab(NULL) +
         xlab(NULL) +
         coords.obj

    color.aes <- .check_aes_exits(p$mapping, c('color', 'colour')) 
    if (!is.null(color.aes)) {
        type.color.value <- p$data |> dplyr::pull(!!color.aes) 
        if (inherits(type.color.value, 'numeric')) {
            p <- p + scale_color_gradient2(low = '#3A3A98', high = '#832424')
        }else{
            p <- p + scale_color_discrete()
        }
    }

    return(p)
})


#' @export
setMethod('plot_point_features', 'SpatialExperiment', 
          function(
            x,
            assayName,
            features = NULL,
            reducedDim.type=NULL,
            dims = c(1, 2),
            mapping = NULL,
            ncol = 3,
            image.plot = TRUE,
            image.rotate.degree = NULL,
            image.mirror.axis = 'h',            
            remove.point = FALSE,
            ...){
    callNextMethod()
})
