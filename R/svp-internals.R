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
.runMCA.internal <- function(x, reduction.name = 'MCA', ncomponents, coords = NULL){
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
  if (!is.null(coords)){
      mca <- cbind(mca, coords)
  }
  colnames(MCA$featuresCoordinates) <- component
  rownames(MCA$featuresCoordinates) <- features.nm
  attr(mca, 'genesCoordinates') <- MCA$featuresCoordinates
  attr(mca, 'stdev') <- SVD$d[-1]
  return(mca) 
}

.build.new.reduced <- function(rd.res, cells, features = NULL){
  
  if (!is.null(features)){
    rd.f.res <- attr(rd.res, 'genesCoordinates')
    rd.f.res <- rd.f.res[features,,drop=FALSE]
  }else{
    rd.f.res <- NULL
  }
  
  rd.res <- rd.res[cells,,drop=FALSE]
  attr(rd.res, "genesCoordinates") <- rd.f.res
  return(rd.res)

}

#' @importFrom BiocNeighbors findKmknn
#' @importFrom BiocParallel SerialParam
.build.knn.graph <- function(x, 
                       knn.k.use = 350, 
                       fun.nm = "findKmknn", 
                       BPPARAM = SerialParam(), 
                       weighted.distance = TRUE,
                       graph.directed = FALSE,
		       ...){
  dots.params <- list(...)
  all.params <- .extract_dot_args(fun.nm, 
				  c('X', 'k', 'BPPARAM'), 
				  dots.params)
  all.params$X <- x
  all.params$k <- knn.k.use
  all.params$BPPARAM <- BPPARAM
  knn.res <- do.call(fun.nm, all.params)
  knn.edges <- .extract_edge(knn.res, x, weighted.distance = weighted.distance)
  knn.graph <- .build.graph(knn.edges, graph.directed = graph.directed) 
  return(knn.graph)
}


.subset_ind <- function(da, nm){
  if (!is.null(nm)){
     x <- intersect(nm, rownames(da))
  }else{
     x <- rownames(da)
  }
  return(x)
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
  knn.edges <- data.frame(knn.edges)
  colnames(knn.edges) <- c('from', 'to')

  if (weighted.distance){
    knn.edges$weight <- .normalize_dist(knn$distance)
  }
  return(knn.edges)
}

#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph simplify
.build.graph <- function(
		  edges, 
		  graph.directed = FALSE){
  g <- igraph::graph_from_data_frame(
         edges, 
	 directed = graph.directed
       ) 
  g <- igraph::simplify(g)
  return(g)
}


.generate.gset.seed <- function(g, gset.idx.list, min.sz = 1, max.sz = Inf, verbose = FALSE){
  
  total.len <- lapply(gset.idx.list, length)

  if (any(total.len <=1) && verbose){
      cli::cli_warn(c("Some gene sets have size one.", 
		      "You've supplied 'min.sz = {min.sz},'
		      consider setting 'min.sz > 1'."))  
  }
  #if (verbose){
  #     cli::cli_inform(c("The gene sets which length lower than {.cls {class(min.sz)}} {min.sz}
  #                       and larger than {.cls {class(max.sz)}} {max.sz} are removed"))
  #}
  gset.idx.list <- gset.idx.list[total.len >= min.sz & total.len < max.sz]

  x <- matrix(0, 
	      nrow = igraph::vcount(g), 
	      ncol = length(gset.idx.list)
       )
  nm <- igraph::vertex_attr(g, 'name')
  rownames(x) <- nm
  if (is.null(names(gset.idx.list))){
     cli::cli_abort('The gene set list must have names.')
  }
  colnames(x) <- names(gset.idx.list)
  
  for (i in seq(length(gset.idx.list))){
     x[,i] <- ifelse(nm %in% gset.idx.list[[i]], 1, 0)
  }
  
  x <- Matrix::Matrix(x, sparse = TRUE)
  return(x)
}

.use.nngp <- function(x, coords, n = 10, order = 'AMMD', BPPARAM = SerialParam()){
  rlang::check_installed('BRISC', "for identify svg using NNGP method.")
  coords <- .normalize.coords(coords)
  # calculate ordering of coordinates
  brisc.order <- BRISC::BRISC_order(coords, order = order, verbose = FALSE)
  # calculate nearest neighbors
  brisc.nn <- BRISC::BRISC_neighbor(coords, n.neighbors = n, n_omp = 1, 
			     search.type = "tree", 
			     ordering = brisc.order, verbose = FALSE) 

  brisc.out <- bplapply(seq(nrow(x)), function(i){
    est.out <- BRISC::BRISC_estimation(coords, y = x[i, ], cov.model = "exponential", 
				       ordering = brisc.order, 
				       neighbor = brisc.nn, verbose = FALSE)
    est.out <- c(est.out$Theta, loglik = out_i$log_likelihood)
    return(est.out)
  }, BPPARAM = BPPARAM)

  brisc.out <- do.call('rbind', brisc.out)

}

.add.int.rowdata <- function(sce, getfun, setfun1, setfun2, namestr, val){
  tmp <- getfun(sce) 
  if (length(tmp) == 0 ){
    sce <- setfun1(sce, value=list())
  }
  sce <- setfun2(x=sce, type=namestr, value=val)
  return (sce)
}

#' @importFrom S4Vectors DataFrame List
.extract.features.score <- function(x, gset.nm, features, gset.idx.list){
  keep.gset <- x[rownames(x) %in% gset.nm, ,drop=FALSE]
  keep.gset <- keep.gset[,colnames(keep.gset) %in% features,drop=FALSE]
  rnm <- rownames(keep.gset)
  cnm <- colnames(keep.gset)
  keep.gset.list <- gset.idx.list[names(gset.idx.list) %in% rownames(keep.gset)]
  res <- ExtractFeatureScoreCpp(keep.gset, rnm, cnm, keep.gset.list)
  res <- DataFrame(features.score = List(res))
  rownames(res) <- rnm
  return(res)
}

.identify.svg <- function(x, 
                          coords, 
                          n = 100, 
                          permutation = 100, 
                          p.adjust.method="fdr",
                          BPPARAM = SerialParam()
			  ){

  coords <- .normalize.coords(coords)
  
  if (BPPARAM$workers == 1){
      res <- CalSpatialKldCpp(coords, x, n, permutation)
  }else{
 
      bgkld <- CalBgSpatialKld(coords, n)

      res <- bplapply(seq(nrow(x)), function(i){
                CalSpatialKld(coords, x[i, ], bgkld, n, permutation) 
             }, BPPARAM = BPPARAM)

      res <- do.call('rbind', res)
  }

  rownames(res) <- rownames(x)
  colnames(res) <- c("sp.kld", "boot.sp.kld.mean", "boot.sp.kld.sd", "sp.kld.pvalue")
  sp.kld <- res[,'sp.kld']
  kld.rank <- rep(NA, length(sp.kld))
  kld.rank[!is.na(sp.kld)] <- rank(sp.kld, na.last=NA)
  res <- cbind(res,
               sp.kld.rank = max(kld.rank, na.rm=TRUE) - kld.rank + 1,
               sp.kld.p.adj = p.adjust(res[,"sp.kld.pvalue"], method = p.adjust.method))
  return(res)

}

.normalize.coords <- function(x){
  range_all <- max(apply(x, 2, function(col) diff(range(col))))
  x <- apply(x, 2, function(col) (col - min(col)) / range_all)
}
