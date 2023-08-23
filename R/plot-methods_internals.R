#' @importFrom ggplot2 ggtitle facet_wrap element_text 
#' @import ggplot2
.feature_setting <- function(features, ncol) {
    if (length(features) == 1) {
        res <- list(ggtitle(features),
            theme(plot.title=element_text(size=rel(1.5), face='bold'))
        )
    }else if(missing(features) || is.null(features)){
        res <- theme_bw2()
    }else{
        res <- list(facet_wrap(~.features, ncol=ncol),
            theme_bw2()
        )
    }
    return(res)
}

##' @importFrom ggplot2 %+replace%
theme_bw2 <- function(...) {
    theme_bw() %+replace%
    theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = 'bold'),
        strip.background = element_rect(color = NA)
      ) %+replace%
    theme(...)
}

.adjust_coords_by_image <- function(coords, x, ...){
    coords <- SpatialExperiment::scaleFactors(x)[1] * coords
    return(coords) 
}


.build_image_annotation <- function(x, rotate.degree = NULL, mirror.axis = NULL, ...){
    img <- getImg(x)
    if (!is.null(rotate.degree)){
        img <- SpatialExperiment::rotateImg(img, rotate.degree)
    }

    if (!is.null(mirror.axis)){
        img <- SpatialExperiment::mirrorImg(img, mirror.axis)
    }
    
    annotation_custom(grob = grid::rasterGrob(imgRaster(img)), 
                      xmin = 1, 
                      ymin = 1, 
                      xmax = dim(img)[2], 
                      ymax = dim(img)[1]
    )

}


