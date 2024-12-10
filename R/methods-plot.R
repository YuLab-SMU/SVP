#' @title plot_heatmap_globalbv
#'
#' @description
#' visulize the result of global bivariate spatial analysis with heatmap
#'
#' @param globalbv the result of \code{runGLOBALBV} with \code{action = 'get'}.
#' @param moran.t the result of global spatial variable features for one type features
#' default is NULL. or \code{runDetectSVG} and then using \code{svDf} to extract the 
#' result.
#' @param moran.l the result of global spatial variable features for anothor type features
#' default is NULL.
#' @param lisa.l the result of \code{cal_lisa_f1} for one type features
#' @param lisa.t the result of \code{cal_lisa_f1} for another type features.
#' @param max.point.size the max point size for main dotplot, default is 4.5.
#' @param font.size the size of font when the triangle heatmap is displayed, default is 2.5.
#' @param limits.size adjust the limit of point size for main dotplot via \code{limits} of 
#' \code{scale_size_continuous}, default is NULL.
#' @param limits.colour adjust the limit of point colour for main dotplot via \code{limits} of
#' \code{scale_fill_gradient2}, default is NULL.
#' @param dist.method the distance measure to be used for the result of global bivariate spatial.
#' which is to measure the dissimilarity between the features, default is \code{'euclidean'}.
#' @param hclust.method the agglomeration method to be used for the result of global bivariate spatial.
#' which is also to measure the similarity between the features, default is \code{'averate'}.
#' @param threshold numeric the threshold to display the point with the significance level,
#' default is 0.05.
#' @return a ggplot2 or aplot object
#' @importFrom ggtree ggtree td_filter layout_dendrogram
#' @importFrom ggfun theme_blinds
#' @examples
#' data(hpda_spe_cell_dec)
#' gbv.res <- runGLOBALBV(
#'               hpda_spe_cell_dec, features1=rownames(hpda_spe_cell_dec), 
#'               assay.type=1, add.pvalue=TRUE, permutation=NULL, alternative='greater'
#'            )
#' 
#' moran.res <- runDetectSVG(hpda_spe_cell_dec, assay.type=1) |> svDf() 
#' 
#' lisa.res <- runLISA(hpda_spe_cell_dec, features=rownames(hpda_spe_cell_dec), assay.type=1)
#' 
#' lisa.f1 <- cal_lisa_f1(hpda_spe_cell_dec, lisa.res, group.by='cluster_domain') 
#' 
#' plot_heatmap_globalbv(gbv.res, moran.t=moran.res, lisa.t=lisa.f1)
#' @export
plot_heatmap_globalbv <- function(globalbv,
                                  moran.t = NULL,
                                  moran.l = NULL,
                                  lisa.t = NULL,
                                  lisa.l = NULL,
                                  max.point.size = 4.5,
                                  font.size = 2.5,
                                  limits.size = NULL,
                                  limits.colour = NULL,
                                  dist.method = "euclidean",
                                  hclust.method = "average",
				  threshold = 0.05
                                  ){
   p1 <- f1 <- p.lisa.l <- p.lisa.t <- p.moran.l <- p.moran.t <- NULL
   flag.square <- FALSE
   pvalnm <- NULL
   if (length(globalbv) == 2 && inherits(globalbv, 'list')){
     index <- names(globalbv)[[1]]
     if (!is.null(globalbv[[2]])){
       pvalnm <- names(globalbv)[[2]]
     }
     flag.square <- identical(rownames(globalbv[[1]]), colnames(globalbv[[1]]))
     globalbv |> as_tbl_df(diag=TRUE, flag.clust=TRUE, dist.method = dist.method, hclust.method=hclust.method) -> da
     p1 <- dist(globalbv[[1]], dist.method) |> 
               hclust(hclust.method) |> 
               ggtree(branch.length='none', ladderize=FALSE)
     if (!flag.square){
        f1 <- p1 + layout_dendrogram()
        p1 <- dist(globalbv[[1]] |> t(), dist.method) |> hclust(hclust.method) |> ggtree(branch.length="none", ladderize=FALSE)
     }
   }else{
     da <- globalbv
     index <- colnames(da)[[3]]
     if (ncol(da)>3){
       pvalnm <- colnames(da)[[4]]
     }
   }
   da <- da |>
         dplyr::mutate(
           group = dplyr::if_else(!!rlang::sym(index) >0, "P", "N")
         )
   p2 <- .plot_main_dot(da, index, pvalnm, max.point.size, limits.size, limits.colour, threshold)
   
   if (!is.null(lisa.l)){
      p.lisa.l <- .plot_lisa_heatmap(lisa.l, coord.flip=TRUE) 
      p2 <- p2 |> aplot::insert_left(p.lisa.l, width = .1)
   }

   if (!is.null(lisa.t)){
      p.lisa.t <- .plot_lisa_heatmap(lisa.t)   
      p2 <- p2 |> aplot::insert_top(p.lisa.t, height = .1)
   }

   if (!is.null(moran.l)){
      p.moran.l <- .plot_moran_barplot(moran.l, coord.flip=TRUE)
      p2 <- p2 |> aplot::insert_left(p.moran.l, width = .18)
   }

   if (!is.null(moran.t)){
      p.moran.t <- .plot_moran_barplot(moran.t)
      p2 <- p2 |> aplot::insert_top(p.moran.t, height = .18)
   }
   if (!is.null(p1)){
       if (flag.square){
           dt <- data.frame(x=levels(da$x), y=levels(da$y))
           layer.text <- geom_text(dt, mapping=aes(x=!!sym("x"),y=!!sym("y"), label=!!sym("y")), 
                                        vjust=0.5, hjust=0, color='black', inherit.aes=FALSE, angle=0, 
                                        size = font.size, nudge_x= 0.5, nudge_y= 0)
           
           if (inherits(p2, 'aplot')){
               p2[[1]] <- p2[[1]] + layer.text + guides(x=guide_axis(angle=45, position='none'), y=guide_axis(position='none')) 
           }else{
               p2 <- p2 + layer.text + guides(x=guide_axis(angle=45, position='none'), y=guide_axis(position='none'))
           }
       }else{
           guide.and.theme <- list(guides(y = guide_axis(position = 'right')), 
                                   theme_blinds(colour=c("gray88", "white")))
           if (inherits(p2, "aplot")){
               p2[[1]] <- p2[[1]] + guide.and.theme
           }else{
               p2 <- p2 + guide.and.theme
           }
       }
       p2 <- p2 |> aplot::insert_left(p1, width=.12)
       if (!is.null(f1)){
           p2 <- p2 |> aplot::insert_top(f1, height = .12)
       }
   }else{
       if (inherits(p2, 'aplot')){
           p2[[1]] <- p2[[1]] + theme_blinds(colour=c("gray88", "white"))
       }else{
           p2 <- p2 + theme_blinds(colour = c("gray88", "white"))
       }
   }   

   return(p2)
}

#' @importFrom ggstar geom_star scale_starshape_manual
#' @importFrom ggplot2 scale_fill_gradient2 scale_size_continuous margin unit 
#' @importFrom ggplot2 aes coord_cartesian guides xlab ylab guide_colorbar
.plot_main_dot <- function(da, 
                           index, 
                           pvalnm,
                           max.point.size, 
                           limits.size, 
                           limits.colour, 
                           threshold
   ){
   ggstar <- "ggstar"
   require(ggstar, quietly=TRUE, character.only=TRUE)
   p <- ggplot(da,
                mapping = aes(x = !!rlang::sym("x"), 
                    y = !!rlang::sym("y"), 
                    fill = !!rlang::sym(index), 
                    size = abs(!!rlang::sym(index)))
                ) +
         geom_star(mapping = aes(starshape = !!rlang::sym("group")), starstroke=.1) +
         scale_starshape_manual(values = c(P=15, N=13), guide='none') +
         scale_fill_gradient2(
           low = scales::muted("blue"),
           high = scales::muted("red"),
           limits = limits.colour,
           guide = guide_colorbar(
                theme = theme(
                  legend.title = element_text(size=7),
                  legend.text = element_text(size=5),
                  legend.key.width  = grid::unit(.5, "lines"),
                  legend.key.height= grid::unit(3, "lines"))
           )
         ) +
         scale_size_continuous(
           range = c(.1, max.point.size),
           limits = limits.size,
           guide = guide_legend(
                override.aes = list(starshape = 15),
                theme = theme(
                  legend.title = element_text(size=7),
                  legend.text = element_text(size=5),
                  legend.key.width  = grid::unit(.4, "cm"),
                  legend.key.height = grid::unit(.4, "cm")
                )
           )
         ) +
         xlab(NULL) +
         ylab(NULL) +
         guides(x=guide_axis(angle=-45)) +
         coord_cartesian(clip= 'off') +
         theme_bw() +
         theme(
            panel.grid = element_blank(),
            panel.border = element_blank(),
            plot.margin = margin(0, 0, 3, 0),
         )
   if (!is.null(pvalnm)){
     p <- p +
          geom_star(
            data = td_filter(!!rlang::sym(pvalnm) <= threshold),
            starshape = 29,
            starstroke =.1
          )           
   }

   return(p)
}


#' @importFrom rlang sym
#' @importFrom dplyr case_when
.tidy_moran_data_for_barplot <- function(da){
    da <- da |>
          tibble::rownames_to_column(var= 'type') |>
          dplyr::mutate(annot = dplyr::case_when(
                            !!sym("padj") <= 0.001 ~ "***",
                            !!sym("padj") > 0.001 & !!sym("padj") <= 0.01 ~ "**",
                            !!sym("padj") > 0.01 & !!sym("padj") <= 0.05 ~ "*",
                            TRUE ~ NA
                         )
          )
    return(da)  
}

#' @importFrom ggplot2 ggplot geom_text geom_col geom_tile theme theme_bw expansion
#' @importFrom ggplot2 scale_y_continuous scale_x_reverse 
.plot_moran_barplot <- function(da, width=4/5, coord.flip = FALSE){
    index <- colnames(da)[[1]]
    da <- .tidy_moran_data_for_barplot(da)
    if (coord.flip){
        mapping <- aes(y = !!sym("type"), x = !!sym(index), fill = !!sym("type"))
    }else{
        mapping <- aes(x = !!sym("type"), y = !!sym(index), fill = !!sym("type"))
    }
    p <- ggplot(da, mapping)
    if (!is.null(width)){
        p <- p + geom_col(show.legend = FALSE, width = width)
    }else{
        p <- p +
            geom_col(show.legend=FALSE)
    }
    p <- p +
            theme_bw() +
            theme(panel.grid = element_blank(),
              panel.border = element_blank(),
           )
    if (coord.flip){
        p <- p +
             geom_text(aes(label = !!sym("annot")), angle = 90, vjust = 1, hjust=0.5) +
             xlab(NULL) +
             ylab(NULL) +
             guides(x = guide_axis(position='bottom'), y=guide_axis(position='none')) +
             theme(axis.line.x.bottom = element_line(), axis.text.x.bottom = element_text())
        p <- p + scale_x_reverse(expand = expansion(mult = c(0.05, 0)))
    }else{
        p <- p +
             geom_text(aes(label = !!sym("annot"))) +
             xlab(NULL) +
             ylab(NULL) +
             guides(x=guide_axis(position='none'), y=guide_axis(position='right')) +
             theme(axis.line.y.right = element_line(), axis.text.y.right = element_text())
        p <- p + scale_y_continuous(expand = expansion(mult = c(0, .05)))
    }
    
    return(p)   
}

.tidy_lisa_f1 <- function(da){
  da |> 
    tibble::as_tibble(rownames='features') |> 
    tidyr::pivot_longer(cols=!"features", names_to="cluster", values_to="F1")

}

#' @importFrom ggplot2 guide_axis element_blank element_text element_line guide_legend
.plot_lisa_heatmap <- function(da, coord.flip=FALSE){
    da <- .tidy_lisa_f1(da)
    if (coord.flip){
        mapping <- aes(y=!!sym("features"), x = !!sym("cluster"), fill = !!sym("cluster"), alpha=!!sym("F1"))
    }else{
        mapping <- aes(x=!!sym("features"), y = !!sym("cluster"), fill = !!sym("cluster"), alpha=!!sym("F1"))
    }
    p <- ggplot(da, mapping) +
         geom_tile(data=td_filter(!!sym("F1")>0), show.legend=c(fill=FALSE, alpha=TRUE))
    p <- p + guides(alpha = guide_legend(theme = theme(legend.key.width = unit(.3, 'cm'),
                                                       legend.key.height = unit(.3, 'cm')))
    )
    if (coord.flip){
        p <- p + guides(y = 'none', x = guide_axis(angle=-50))
    }else{
        p <- p + guides(x = 'none', y = guide_axis(position = 'right'))
    }
    p <- p +
         theme_bw() +
         theme(panel.grid = element_blank()) +
         xlab(NULL) +
         ylab(NULL)
    return(p)
}
