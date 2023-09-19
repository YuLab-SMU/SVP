#' @importFrom BiocParallel bplapply MulticoreParam
.cal_cor <- function(x, img, coords, beta = 2, 
                     cor.method = 'pearson', 
                     BPPARAM = NULL, threads = 5, p.adjust.method = 'fdr', 
                     index.image = 1, k = 999, ...){
    if (is.null(BPPARAM)) {
        BPPARAM <- MulticoreParam(workers = threads)
    }
    if (nrow(img) == 0){
        return(NULL)
    }else{
        if (is.numeric(index.image)){
            index.image <- unique(img$image_id)[index.image]
        }
        img <- img[img$image_id == index.image,]
        coords <- coords * img[['scaleFactor']]
        img <- as.raster(img[['data']][[1]])
    }
    coords.dist <- as.matrix(dist(coords))
    color.features <- extract_color(img, beta = beta, coords)
    res <- bplapply(seq_len(nrow(x)), function(i){
        xi <- x[i,]
        tmp <- cor.test(as.numeric(xi), as.numeric(color.features), method = cor.method)
        #tmp <- crossCorrelation(as.numeric(xi), as.numeric(color.features), w = coords.dist, k = k)
        #res.i <- c(pvalue.cor.img = tmp$global.p, value.cor.img = as.numeric(tmp$I))
	res.i <- c(pvalue.cor.img = tmp$p.value, value.cor.img = as.numeric(tmp$estimate))
        return(res.i)
        }, 
        BPPARAM = BPPARAM
    )
    res <- do.call('rbind', res)
    p.adj.cor.img <- p.adjust(res[,1], method = p.adjust.method)
    res <- cbind(res, p.adj.cor.img)
    rownames(res) <- rownames(x)
    return(res)
}

## This function is from CelliD, because it has many depends.
## to better developement, it was coped and modified diretcly.
.runMCA.internal <- function(x, reduction.name = 'MCA', ncomponents){
  rlang::check_installed('irlba', '`runMCA()`')
  if (inherits(x, 'dgTMatrix')){
     x <- as(x, "CsparseMatrix")
  }
  x <- x[rowVars(x) != 0, ]
  cells.nm <- colnames(x)
  features.nm <- rownames(x)
  if (ncol(x) < 4){
     cli::cli_abort(c('The length of variables is', ncol(x), ', which is too few to runMCA().'))
  }
  ncomponents <- min(ncol(x) - 2, ncomponents)
  x <- as.matrix(x)
  message("Computing Fuzzy Matrix")
  MCAPrepRes <- MCAStep1(x)
  message("Computing SVD")
  SVD <- irlba::irlba(A = MCAPrepRes$Z, nv = ncomponents + 1, nu = 1)[seq(3)]
  message("Computing Coordinates")
  MCA <- MCAStep2(Z = MCAPrepRes$Z, V = SVD$v[, -1], Dc = MCAPrepRes$Dc)
  component <- paste0(reduction.name, "_", seq(ncol(MCA$cellsCoordinates)))
  mca <- MCA$cellsCoordinates
  colnames(mca) <- component
  rownames(mca) <- cells.nm
  colnames(MCA$featuresCoordinates) <- component
  rownames(MCA$featuresCoordinates) <- features.nm
  attr(mca, 'featuresCoords') <- MCA$featuresCoordinates
  attr(mca, 'stdev') <- SVD$d[-1]
  return(mca) 
}


#' @importFrom BiocNeighbors findKmknn
#' @importFrom BiocParallel SerialParam
.build.knn <- function(x, 
		       knn.k.use = 300, 
		       fun.nm = "findKmknn", 
		       BPPARAM = SerialParam(), 
		       weighted.distance = TRUE,
		       ...){
  dots.params <- list(...)
  all.params <- .extract_dot_args(fun.nm, 
				  c('X', 'k', 'BPPARAM'), 
				  dots.params)
  all.params$X <- x
  all.params$k <- knn.k.use
  all.params$BPPARAM <- BPPARAM
  knn.res <- do.call(fun.nm, all.params)
  knn.edges <- .extract_edge(knn.res, weighted.distance = weighted.distance)
  

}

.subset_data <- function(x, n){
  if (!is.null(n)){
    x <- x[n,,drop=FALSE]
  }
  return(x)
}

.extract_dot_args <- function(fun.nm, used.arg, dots){
  total.args <- names(as.list(args(fun.nm))) 
  left.args <- setdiff(c(total.args, "distance"), c(used.arg, "...", ""))
  dots <- dots[intersect(names(dots), left.args)]
  return(dots)  
}

.normalize_dist <- function(x){
  x <- as.vector(x)
  (x - min(x))/(max(x) - min(x))
}

.extract_edge <- function(knn, x, weighted.distance = TRUE){
  nm <- rownames(x)
  knn.index <- knn$index
  knn.index <- cbind(rep(seq(nrow(knn.index)), 
			 times=ncol(knn.index)), 
		     as.vector(knn.index))
  knn.edges <- cbind(nm[knn.index[,1]], nm[knn.index[,2]])
  colnames(knn.edges) <- c('from', 'to')
  if (weighted.distance){
    knn.edges <- cbind(knn.edges, .normalize_dist(knn$distance))
    colnames(knn.edges)[3] <- 'weight'
  }
  return(knn.edges)
}


.build.graph <- function(edges, 
			 graph.directed = FALSE,
			 weighted.distance = TRUE){
  g <- igraph::graph_from_data_frame(edges, directed = graph.directed) |>
       igraph::simplify()
  return(g)
}
