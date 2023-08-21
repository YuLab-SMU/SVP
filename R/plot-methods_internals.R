#' @importFrom ggplot2 ggtitle facet_wrap element_text 
#' @import ggplot2
.feature_setting <- function(features, ncol) {
    if (length(features) == 1) {
        res <- list(ggtitle(features),
            theme(plot.title=element_text(size=rel(1.5), face='bold'))
        )
    } else {
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
