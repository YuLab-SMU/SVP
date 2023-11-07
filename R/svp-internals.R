## This function is from CelliD, because it has many depends.
## to better developement, it was coped and modified diretcly.
.runMCA.internal <- function(x, reduction.name = 'MCA', ncomponents, coords = NULL){
  rlang::check_installed(c('irlba', 'Matrix', 'matrixStats'), '`runMCA()`')
  if (inherits(x, 'dgTMatrix')){
     x <- as(x, "CsparseMatrix")
  }else if (inherits(x, 'matrix')){
     x <- Matrix::Matrix(x, sparse = TRUE)
  }
  #x <- x[matrixStats::rowVars(x) != 0, ]
  cells.nm <- colnames(x)
  features.nm <- rownames(x)
  if (ncol(x) < 4){
     cli::cli_abort(c('The length of variables is', ncol(x), ', which is too few to runMCA().'))
  }
  ncomponents <- min(ncol(x) - 2, ncomponents)
  #x <- as.matrix(x)
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
  attr(mca, 'mcaFlag') <- TRUE
  attr(mca, 'stdev') <- SVD$d[-1]
  return(mca) 
}

.build.new.reduced <- function(rd.res, cells, features = NULL, rd.f.nm = 'genesCoordinates'){
  rd.old <- rd.res
  if (!is.null(features)){
    rd.f.res <- attr(rd.res, rd.f.nm)
    rd.f.res <- rd.f.res[features,,drop=FALSE]
  }else{
    rd.f.res <- NULL
  }
  
  rd.res <- rd.res[cells,,drop=FALSE]
  attr(rd.res, rd.f.nm) <- rd.f.res

  othernm <- setdiff(names(attributes(rd.old)), names(attributes(rd.res)))
  if (length(othernm) > 0){
      for (i in seq(length(othernm))){
          attr(rd.res, othernm[i]) <- attr(rd.old, othernm[i])
      }
  }

  return(rd.res)
}

#' @importFrom BiocNeighbors findKmknn
#' @importFrom BiocParallel SerialParam
.build.knn.graph <- function(x, 
                       knn.k.use = 350, 
                       fun.nm = findKmknn, 
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


.filter.gset.gene <- function(x, gset.idx.list, min.sz = 10, max.sz = Inf, 
			      gene.occurrence.rate = .4){
    
    total.genes <- x
    gset.gene.num <- lapply(gset.idx.list, function(i)length(unique(i))) |> unlist()
    exp.gene.num  <- lapply(gset.idx.list, function(i) sum(unique(i) %in% total.genes)) |> unlist()
    if (any(gset.gene.num <=1) && min.sz == 1){
        cli::cli_warn(c("Some gene sets have size one.",
                        "You've supplied 'min.sz = {min.sz},'
                        consider setting 'min.sz > 1'."))
    }    
    
    gene.num <- data.frame(
                  exp.gene.num = exp.gene.num, 
	          gset.gene.num = gset.gene.num, 
	          gene.occurrence.rate = exp.gene.num/gset.gene.num
                )
    rownames(gene.num) <- names(gset.idx.list)
    gene.num <- gene.num[gene.num$gset.gene.num >= min.sz & gene.num$gene.occurrence.rate >= gene.occurrence.rate,]

    return(as.matrix(gene.num))

}


.generate.gset.seed <- function(g, 
                                gset.idx.list
				){
  
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

#' @importFrom rlang .data
#' @importFrom stats p.adjust
#' @importFrom BiocParallel bplapply
.identify.svg <- function(x, 
                          coords, 
                          n = 100, 
                          permutation = 100, 
                          p.adjust.method="fdr",
                          BPPARAM = SerialParam(),
			  random.seed = 123
			  ){

  rlang::check_installed("withr", "is required to reproducible in the identification of SVG or SVP.")

  coords <- .normalize.coords(coords)
  
  lims <- c(range(coords[,1]), range(coords[,2]))
  h <- c(ks::hpi(coords[,1]), ks::hpi(coords[,2]))
 
  gx <- seq.int(lims[1], lims[2], length.out = n)
  gy <- seq.int(lims[3], lims[4], length.out = n)

  indx <- findIntervalCpp(coords[, 1], gx)
  indy <- findIntervalCpp(coords[, 2], gy)

  axm <- outergrid(gx, coords[, 1]) 
  aym <- outergrid(gy, coords[, 2])
 
  bgkld <- CalBgSpatialKld(coords, axm, aym, h, indx, indy)

  res <- bplapply(seq(nrow(x)), function(i){
            withr::with_seed(random.seed, CalSpatialKld(x[i, ], bgkld, axm, aym, h, indx, indy, permutation, random.seed))
         }, BPPARAM = BPPARAM)

  res <- do.call('rbind', res)

  rownames(res) <- rownames(x)
  colnames(res) <- c("sp.kld", "boot.sp.kld.mean", "boot.sp.kld.sd", "sp.kld.pvalue")
  res <- cbind(res,
	       sp.kld.p.adj = p.adjust(res[,"sp.kld.pvalue"], method = p.adjust.method)
            ) |> data.frame()
  res <- res |> dplyr::arrange(.data$sp.kld.p.adj, dplyr::desc(.data$sp.kld)) |>
         dplyr::mutate(sp.kld.rank = seq(nrow(res))) |> as.matrix()
  return(res)
}

.normalize.coords <- function(x){
  range_all <- max(apply(x, 2, function(col) diff(range(col))))
  x <- apply(x, 2, function(col) (col - min(col)) / range_all)
}
