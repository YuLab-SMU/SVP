## This function is from CelliD, because it has many depends.
## to better developement, it was coped and modified diretcly.
.runMCA.internal <- function(x, reduction.name = 'MCA', ncomponents, coords = NULL){
  rlang::check_installed(c('irlba', 'Matrix'), '`runMCA()`')
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

pairDist <- function(x, y){
    z <- fastPDist(y, x)
    z[is.na(z)] <- 6.629066e-15
    rownames(z) <- rownames(y)
    colnames(z) <- rownames(x)
    return(z)
}

# #' build the adjacency matrix with the distance matrix
# #' @param x a dist matrix
# #' @param top.n number of the nearest neighbors.
# #' @param weighted.distance whether consider the distance between nodes.
# #' @return adjacency matrix

#' @importFrom Matrix sparseMatrix Matrix
.build.adj.m <- function(x, top.n = 600, weighted.distance = FALSE){
    rnm <- rownames(x)
    cnm <- colnames(x)
    
    top.n <- min(nrow(x), top.n)
    i <- apply(x, 2, order)[seq(top.n), ]
    j <- rep(seq(ncol(x)), each = top.n)
    adj.m <- Matrix::sparseMatrix(i, j, x = 1, dims = c(length(rnm),
        length(cnm)), dimnames = list(rnm, cnm))
    if (weighted.distance){
        adj.m <- adj.m * x
    }
    return(adj.m)
}

#' @importFrom DelayedMatrixStats colRanks
.build.adj.m_by_expr <- function(x, top.n = 600, weighted.distance = FALSE){
    expr.rank.dist <- suppressWarnings(1 - (DelayedMatrixStats::colRanks(x, preserveShape=TRUE, useNames=TRUE) / (nrow(x) + 1)))
    adj.m <- .build.adj.m(expr.rank.dist, top.n, weighted.distance)
    if (weighted.distance){
        adj.m@x <- .normalize_dist(adj.m@x)
    }
    return(adj.m)
}

.join.adj.m <- function(fs2cell, cell2cell, fs2fs){
    x12 <- Matrix::rbind2(cell2cell, fs2cell)
    x13 <- Matrix::cbind2(fs2cell, fs2fs)
    x <- Matrix::cbind2(x12, Matrix::t(x13))
    return(x)
}

.calculate.odds.weighted <- function(x, y, z){
    x2 <- x
    x@x <- rep(1, length(x@x))
    q2 <- Matrix::t(x) %*% (1 - y)
    w1 <- (Matrix::t(x2) %*% (1 - y))/q2
    w2 <- (Matrix::t(x2) %*% y) / as.matrix(z+1)
    res <- as.matrix(w1/w2)
    res[is.na(res)] <- 1e-10
    res[res==0] <- 1e-10
    as.data.frame(res, check.names=FALSE)
}

.internal.pWNCHypergeo <- function(x, m1, m2, n, odds, lower.tail=FALSE){
    mapply(BiasedUrn::pWNCHypergeo, x = x, m1 = m1, m2 = m2, n = n, odds = odds,
           MoreArgs=list(lower.tail=lower.tail))
}

# #' Perform hyper geometric test on cells with the nearest neighbors features
# #' @param fs2cell.adj the adjacency matrix of gene to cell.
# #' @param fs2gset.adj the adjacency matrix of gene to gene set (or pathway).
# #' @param gene.set.list the pathway or gene set list
# #' @param top.n the number of nearest neighbors features
# #' @return a matrix, the result of hypergeometric for each cell to each pathway.

#' @importFrom stats phyper
#' @importFrom fastmatch %fin%
.run_hgt <- function(fs2cell.adj, 
                     fs2gset.adj, 
                     gene.set.list, 
                     cells.gene.num = NULL,
		     m = NULL,
                     top.n = 600, 
                     combined.cell.feature = FALSE,
                     weighted.distance = FALSE,
                     method = c("Wallenius", "Hypergeometric")
                    ){
    method <- match.arg(method)
    if (!weighted.distance){
        method <- 'Hypergeometric'
    }
    fs2cell.adj2 <- fs2cell.adj 
    if (weighted.distance){
        fs2cell.adj@x <- rep(1, length(fs2cell.adj@x))
    }
    q <- as.data.frame(suppressWarnings(as.matrix(Matrix::t(fs2cell.adj) %*% fs2gset.adj)) - 1)

    if (is.null(m)){
        m <- vapply(gene.set.list, function(x) sum(x %fin% rownames(fs2cell.adj)), numeric(1))
    }
    
    if (is.null(cells.gene.num)){
        n <- nrow(fs2cell.adj) - m
    }else{
        n <- vapply(cells.gene.num, function(i)i-m, numeric(length(m))) |> t() |> as.data.frame(check.names=FALSE)
    }
    k <- top.n
    if (combined.cell.feature){
        k <- Matrix::colSums(fs2cell.adj)
    }
    if (method == 'Hypergeometric'){
        res <- mapply(phyper, q = q, m = m, n = n, MoreArgs=list(k=k, lower.tail = FALSE))
    }else if (method == 'Wallenius'){
        rlang::check_installed('BiasedUrn', 'for `Wallenius` method.')
        odds <- .calculate.odds.weighted(fs2cell.adj2, fs2gset.adj, q)
        res <- mapply(.internal.pWNCHypergeo, x = q, m1 = m, m2 = n, n = k,
                      odds = odds, MoreArgs=list(lower.tail = FALSE)) 
    }
    rownames(res) <- rownames(q)
    colnames(res) <- colnames(fs2gset.adj)
    res[res==0] <- 1e-10
    res <- -log10(t(res))
    res[res < 0 ] <- 0
    #res <- Matrix::Matrix(res, sparse = TRUE)
    return(res)
}

.build.nndist.graph <- function(
                            cells.rd, 
                            features.rd,
                            cells.sp.coord = NULL,
                            knn.consider.spcoord = FALSE,
                            sp.alpha.add.weight = .2,
                            sp.beta.add.mp.weight = .1,
                            top.n = 600,
                            combined.cell.feature = FALSE, 
                            weighted.distance = FALSE,
                            graph.directed = FALSE,
                            normalize.dist = TRUE,
                            BPPARAM = BiocParallel::MulticoreParam(workers = 3)
                       ){
    if (!combined.cell.feature){
        # This is split the cells and features to build knn
        x <- pairDist(cells.rd, features.rd)
        
        #cell.dist <- pairDist(cells.rd, cells.rd)
        
        cell.dist <- .fusion.sc.sp.dist(
                       cells.rd, 
                       cells.sp.coord, 
                       knn.consider.spcoord, 
                       sp.alpha.add.weight, 
                       sp.beta.add.mp.weight
        )

        fs.dist <- pairDist(features.rd, features.rd)
        top.n <- min(top.n, nrow(features.rd))
        top.n.cell <- min(max(50, round(top.n/10)), nrow(cells.rd)) + 1
        top.n.fs <- min(max(50, round(top.n/10)), nrow(features.rd)) + 1
        adj.m.list <- BiocParallel::bpmapply(.build.adj.m, list(x, cell.dist, fs.dist), 
                                             list(top.n, top.n.cell, top.n.fs), 
                                             MoreArgs=list(weighted.distance),
                                             BPPARAM = BPPARAM
        ) 
        adj.m <- do.call(.join.adj.m, adj.m.list)
        Matrix::diag(adj.m) <- 0
    }else{
        # build knn by merge the MCA space of cells and features 
        total.rd <- rbind(cells.rd, features.rd)
        #total.dist <- pairDist(total.rd, total.rd)
        top.n <- min(top.n, nrow(total.rd)) 
        knn.graph <- .build.knn.graph(total.rd, top.n, 
                                      fun.nm = BiocNeighbors::findKmknn, 
                                      weighted.distance = weighted.distance)
        adj.m <- .extract.adj.m(knn.graph, edge.attr = 'weight') 
        #adj.m <- .build.adj.m(total.dist, top.n, weighted.distance)
    }
    if (weighted.distance){
       if (normalize.dist){
          adj.m@x <- .normalize_dist(adj.m@x)
       }
    }
    return(adj.m)
}

.fusion.sc.sp.dist <- function(cells.rd, 
                               cells.sp.coord, 
                               consider.spcoord = FALSE,
                               alpha = 0.2, 
                               beta=0.1){
    cell.dist <- pairDist(cells.rd, cells.rd)
    if (consider.spcoord && !is.null(cells.sp.coord)){
        cli::cli_inform('Considering the spatial location ...')
        cell.sp.dist <- pairDist(cells.sp.coord, cells.sp.coord)
        nm <- rownames(cell.dist)
        cell.dist <- fusiondist(cell.dist, cell.sp.dist, alpha, beta)
        rownames(cell.dist) <- nm
        colnames(cell.dist) <- nm
    }
    return(cell.dist)
}

#' @importFrom BiocParallel SerialParam
.build.knn.graph <- function(x, 
                       knn.k.use = 600, 
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
  knn.res <- suppressWarnings(do.call(fun.nm, all.params))
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

.normalize.single.score <- function(x){
    if (inherits(x, 'matrix')){
        y <- .normalize_dist(x, beta=0, reverse=FALSE)
        attr(y, 'dim') <- attr(x, 'dim')
        attr(y, 'dimnames') <- attr(x, 'dimnames')
        return(y)
    }else if (inherits(x, 'dgCMatrix')){
        x@x <- .normalize_dist(x@x, beta=0, reverse=FALSE)
    }
    return(x)
}

.normalize_dist <- function(x, beta=.01, reverse=TRUE){
  x <- as.vector(x)
  if (stats::var(x, na.rm = TRUE) == 0){
     return(x)
  }
  x <- beta + (1 - 2 * beta) * (x - min(x))/(max(x) - min(x))
  if (reverse){
     x <- 1 - x
  }
  return(x)
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
    exp.gene.num  <- lapply(gset.idx.list, function(i) sum(unique(i) %fin% total.genes)) |> unlist()
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
  if (inherits(g, 'igraph')){
      num.row <- igraph::vcount(g)
      nm <- igraph::vertex_attr(g, 'name')
  }else{
      num.row <- nrow(g)
      nm <- rownames(g)
  }
  
  if (is.null(names(gset.idx.list))){
     cli::cli_abort('The gene set list must have names.')
  }
  x <- vapply(gset.idx.list, function(i) ifelse(nm %fin% i, 1, 0), numeric(num.row))
  rownames(x) <- nm

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
.extract.features.rank <- function(x, y, features, gset.idx.list){
  y <- y[features %in% rownames(y), ,drop=FALSE]
  keep.gset <- corCpp(Matrix::t(x), Matrix::t(y))
  rownames(keep.gset) <- rownames(x)
  colnames(keep.gset) <- rownames(y)
  keep.gset.list <- gset.idx.list[names(gset.idx.list) %in% rownames(keep.gset)]
  res <- ExtractFeatureScoreCpp(keep.gset, rownames(keep.gset), colnames(keep.gset), keep.gset.list)
  res <- DataFrame(features.score = List(res))
  rownames(res) <- rownames(keep.gset)
  return(res)
}

.identify.svg.by.autocorrelation <- function(
  x, 
  coords, 
  method = 'moransi',
  permutation = 1, 
  scaled = FALSE,
  alternative = 'two.sided',
  p.adjust.method = 'BH',
  ...
  ){
  
  coords <- .normalize.coords(coords)
  coords.dist <- pairDist(coords, coords)
  
  if (method == 'moransi'){
      if (is.null(permutation)){permutation <- 1}
      res <- CalMoransiParallel(x, coords.dist, scaled, permutation)
      colnames(res) <- c('obs', 'expect.moransi', 'sd.moransi', 'pvalue')
  }else if (method == 'gearysc'){
      if (is.null(permutation)){permutation <- 100}
      res <- CalGearyscParallel(x, coords.dist, permutation)
      colnames(res) <- c('obs', 'expect.gearysc', 'sd.gearysc', 'pvalue')
  }

  if (alternative == "two.sided")
      res[,4] <- ifelse(res[,1] <= res[,2], 2 * res[,4], 2 * (1 - res[,4]))
  if (alternative == "greater")
      res[,4] <- 1 - res[,4]

  rownames(res) <- rownames(x)
  res <- cbind(res,
               padj = p.adjust(res[, "pvalue"], method = p.adjust.method)
            ) |> as.data.frame(check.names=FALSE)
  res <- res |> dplyr::arrange(.data$padj, .data$obs) |>
         dplyr::mutate(rank = seq(nrow(res))) |> as.matrix()
  res <- res[match(rownames(x), rownames(res)), ]  
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
  colnames(res) <- c("sp.kld", "boot.sp.kld.mean", "boot.sp.kld.sd", "pvalue")
  res <- cbind(res,
               padj = p.adjust(res[, "pvalue"], method = p.adjust.method)
            ) |> data.frame()
  res <- res |> dplyr::arrange(.data$padj, dplyr::desc(.data$sp.kld)) |>
         dplyr::mutate(rank = seq(nrow(res))) |> as.matrix()
  res <- res[match(rownames(x), rownames(res)), ]
  return(res)
}

.normalize.coords <- function(x){
  range_all <- max(apply(x, 2, function(col) diff(range(col))))
  x <- apply(x, 2, function(col) (col - min(col)) / range_all)
}
