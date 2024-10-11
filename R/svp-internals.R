## This function refers to CelliD, However, we utilize Spectra and Eigen to expedite the 
## computation of SVD and the calculation of feature coordinates (1.5 faster than RunMCA of CelliD).
.runMCA.internal <- function(x, y = NULL, reduction.name = 'MCA', ncomponents, coords = NULL, ...){
  rlang::check_installed(c('RSpectra', 'Matrix'), 'for `runMCA()`')
  if (inherits(x, 'dgTMatrix')){
     x <- as(x, "CsparseMatrix")
  }else if (inherits(x, 'matrix')){
     x <- Matrix::Matrix(x, sparse = TRUE)
  }
  x <- x[DelayedMatrixStats::rowVars(x) != 0, ]
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
  SVD <- RSpectra::svds(MCAPrepRes$Z, k=ncomponents + 1, nu = 1, opts = list(tol = 1e-5))
  V <- SVD$v[,-1,drop=FALSE]
  if (!is.null(y)){
      rlang::check_installed("harmony", 'for `runMCA()` when `group.by.vars` is not NULL.')
      cli::cli_inform("Running Harmony to remove the batch effect")
      group.vars <- unique(colnames(y))
      V <- harmony::RunHarmony(V, y, group.vars, verbose=FALSE, ...)
  }
  message("Computing Coordinates")
  MCA <- MCAStep2(Z = MCAPrepRes$Z, V = V, Dc = MCAPrepRes$Dc)
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
    i <- ParallelColOrder(x, top.n)
    j <- rep(seq(ncol(x)), each = top.n)
    adj.m <- Matrix::sparseMatrix(i, j, x = 1, dims = c(length(rnm),
        length(cnm)), dimnames = list(rnm, cnm))
    if (weighted.distance){
        adj.m <- SpMatElemMultiMat(adj.m, x)
        rownames(adj.m) <- rnm
        colnames(adj.m) <- cnm
    }
    return(adj.m)
}

#' @importFrom DelayedMatrixStats colRanks
.build.adj.m_by_expr <- function(x, top.n = 600, weighted.distance = FALSE){
    expr.rank.dist <- 1 - (DelayedMatrixStats::colRanks(x, useNames=FALSE, preserveShape=TRUE) / (nrow(x) + 1))
    rownames(expr.rank.dist) <- rownames(x)
    colnames(expr.rank.dist) <- colnames(x)
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
    tic()
    cli::cli_inform("ORA analysis ...")
    if (method == 'Hypergeometric'){
        res <- mapply(phyper, q = q, m = m, n = n, MoreArgs=list(k=k, lower.tail = FALSE))
    }else if (method == 'Wallenius'){
        rlang::check_installed('BiasedUrn', 'for `Wallenius` method.')
        odds <- .calculate.odds.weighted(fs2cell.adj2, fs2gset.adj, q)
        res <- mapply(.internal.pWNCHypergeo, x = q, m1 = m, m2 = n, n = k,
                      odds = odds, MoreArgs=list(lower.tail = FALSE)) 
    }
    toc()
    res[res==0] <- 1e-10
    res <- -log10(t(res))
    res[res < 0 ] <- 0
    rownames(res) <- names(gene.set.list)
    colnames(res) <- rownames(q)
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
        adj.m <- .build.knn.adj(total.rd, top.n, 
                                  fun.nm = findKmknn, 
                                  weighted.distance = weighted.distance
                                  ) 
        #adj.m <- .extract.adj.m(knn.graph, edge.attr = 'weight') 
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
#' @importFrom BiocNeighbors findKmknn
.build.knn.adj <- function(x, 
                       knn.k.use = 600, 
                       fun.nm = findKmknn, 
                       BPPARAM = SerialParam(), 
                       weighted.distance = TRUE,
                       bycol = TRUE,
                       ...){
  dots.params <- list(...)
  all.params <- .extract_dot_args(fun.nm, 
                                  c('X', 'k', 'BPPARAM'), 
                                  dots.params)
  all.params$X <- x
  all.params$k <- knn.k.use
  all.params$get.distance <- weighted.distance
  all.params$BPPARAM <- BPPARAM
  knn.res <- suppressWarnings(do.call(fun.nm, all.params))
  res <- .extract_adj(knn.res, rownames(x), weighted.distance, bycol)
  return(res)
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

.extract_adj <- function(knn, x, weighted.distance = TRUE, bycol =TRUE){
  nn <- length(x)
  ind1 <- rep(seq(nn), each=ncol(knn$index))
  ind2 <- c(t(knn$index))
  xx <- ifelse(weighted.distance, c(t(knn$distance)), 1)
  if (bycol){
    res <- Matrix::sparseMatrix(i = ind2, j = ind1, x = xx, dims = c(nn, nn))
  }else{
    res <- Matrix::sparseMatrix(i = ind1, j = ind2, x = xx, dims = c(nn, nn))  
  }
  rownames(res) <- colnames(res) <- x
  return(res) 
}

#.build.graph <- function(
#                  edges, 
#                  graph.directed = FALSE){
#  g <- igraph::graph_from_data_frame(
#         edges, 
#         directed = graph.directed
#       ) 
#  g <- igraph::simplify(g)
#  return(g)
#}


.filter.gset.gene <- function(x, gset.idx.list, min.sz = 10, max.sz = Inf, 
                              gene.occurrence.rate = .4){
    
    total.genes <- x
    gset.gene.num <- lapply(gset.idx.list, function(i)length(unique(i))) |> unlist()
    exp.gene <- lapply(gset.idx.list, function(i)i[unique(i) %fin% total.genes])
    exp.gene.num <- lapply(exp.gene, function(i)length(unique(i))) |> unlist()    
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
    gene.num$geneSets <- exp.gene
    rownames(gene.num) <- names(gset.idx.list)
    gene.num <- gene.num[gene.num$gset.gene.num >= min.sz & gene.num$gene.occurrence.rate >= gene.occurrence.rate,]
    if (nrow(gene.num)==0){
        cli::cli_abort(c("All gene set list was removed since they did not meet these conditions: ",
                         "{.var min.sz} >= {min.sz} and {.var gene.occurrence.rate} >= {gene.occurrence.rate}"), 
                          call = NULL)
    }

    return(gene.num)
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
  ind <- lapply(gset.idx.list, function(i) which(nm %fin% i)) 
  x <- Matrix::sparseMatrix(
         i = ind |> unlist(),
         j = lapply(seq(length(ind)), function(i)rep(i, length(ind[[i]]))) |> unlist(), 
         x = 1,
         dims = c(num.row, length(gset.idx.list))
  )
  rownames(x) <- nm
  colnames(x) <- names(gset.idx.list)
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
  #keep.gset <- corCpp(Matrix::t(x), Matrix::t(y))
  keep.gset <- fast_cor(x, y, method='spearman')
  keep.gset <- keep.gset$r
  #rownames(keep.gset) <- rownames(x)
  #colnames(keep.gset) <- rownames(y)
  keep.gset.list <- gset.idx.list[names(gset.idx.list) %in% rownames(keep.gset)]
  res <- ExtractFeatureScoreCpp(keep.gset, rownames(keep.gset), colnames(keep.gset), keep.gset.list)
  res <- DataFrame(features.score = List(res))
  rownames(res) <- rownames(keep.gset)
  return(res)
}

.identify.svg.by.autocorrelation <- function(
  x,
  coords,
  weight = NULL,
  weight.method = c("voronoi", "knn", "none"), 
  method = 'moransi',
  permutation = 1, 
  scaled = FALSE,
  alternative = c('two.sided', 'greater', 'less'),
  p.adjust.method = 'BH',
  random.seed = 1024,
  ...
  ){

  rlang::check_installed("withr", "is required to be reproducible in the identification of SVG or SVP.")  
  
  params <- list(...)
 
  if ('alternative' %in% names(params)){
      alternative <- params$alternative
      params$alternative <- NULL
  }else{
      if (method == 'gearysc'){
          alternative <- "less"
      }
      if (method %in% c('moransi', "getisord")){
          alternative <- "greater"
      }
  }

  alternative <- match.arg(alternative)
  
  if ('scaled' %in% names(params)){
      scaled <- params$scaled
      params$scaled <- NULL
  }

  if (alternative == 'greater'){
      lower.tail <- 0
  }else{
      lower.tail <- 1
  }

  wm <- .obtain.weight(coords, weight, weight.method, ...)
  if (is.null(permutation)){
      permutation <- 1
  }
  
  if (method == 'moransi'){
      res <- withr::with_seed(random.seed, CalMoransiParallel(x, wm, scaled, permutation, lower.tail))
      colnames(res) <- c('obs', 'expect.moransi', 'sd.moransi', "Z.moransi", 'pvalue')
  }else if (method == 'gearysc'){
      res <- withr::with_seed(random.seed, CalGearyscParallel(x, wm, permutation, lower.tail))
      colnames(res) <- c('obs', 'expect.gearysc', 'sd.gearysc', "Z.gearysc", 'pvalue')
  }else if (method == 'getisord'){
      res <- CalGetisOrdParallel(x, wm, lower.tail)
      colnames(res) <- c("obs", "expect.G", "sd.G", "Z.G", "pvalue")
  }
  
  if (alternative == "two.sided")
      res[,4] <- ifelse(res[,1] <= res[,2], 2 * res[,4], 2 * (1 - res[,4]))

  rownames(res) <- rownames(x)
  res <- cbind(res,
               padj = p.adjust(res[, "pvalue"], method = p.adjust.method)
            ) |> as.data.frame(check.names=FALSE)
  if (method != 'gearysc'){
      res <- res |> dplyr::arrange(.data$padj, .data$pvalue, dplyr::desc(abs(.data$obs))) 
  }else{
      res <- res |> dplyr::arrange(.data$padj, .data$pvalue, dplyr::desc(abs(.data$obs - 1)))
  }
  res <- res |>
         dplyr::mutate(rank = seq(nrow(res)))
  res <- res[match(rownames(x), rownames(res)), ]  
  return(res)
}


.convert_to_distmt <- function(x){
  if (inherits(x, "Graph")){
      res <- .convert_to_distmt.Graph(x)
  }  
  if (inherits(x, 'nb') && !inherits(x, "listw")){
      res <- .convert_to_distmt.nb(x)
  }
  if (inherits(x, "listw")){
      res <- .convert_to_distmt.listw(x)
  }
  if (inherits(x, 'deldir')){
      res <- .convert_to_distmt.deldir(x)
  }
  return(res)
}

.convert_to_distmt.nb <- function(x, w = NULL){
  times <- unlist(lapply(x, length))
  i <- rep(seq(length(x)), times = times)
  j <- unlist(x)
  ind <- j!=0
  j <- j[ind]
  i <- i[ind]
  n <- length(x)
  region_id <- attr(x, "region.id")
  if (is.null(w)){
      w <- rep(1/times, times = times) 
      w <- w[ind]
  }
  res <- sparseMatrix(i = i, j = j, x = w, 
                      dims = rep(n, 2), 
                      dimnames = list(region_id, region_id)) 
  return(res)
}

.convert_to_distmt.listw <- function(x, w = NULL){
  w <- unlist(x$weights)
  res <- .convert_to_distmt.nb(x$neighbours, w = w)
  return(res)
}

.convert_to_distmt.Graph <- function(x){
  w <- unname(table(x$from))
  w <- rep(1/w, times = w)
  region_id <- x$x |> attr("names")
  res <- sparseMatrix(i=x$from, 
                      j=x$to, x=w, 
                      dims=rep(x$np, 2), 
                      dimnames = list(region_id, region_id)) 
  return(res)
}

.convert_to_distmt.deldir <- function(x){
  res <- sparseMatrix(
           i = c(x$delsgs[,5], x$delsgs[,6]), 
           j = c(x$delsgs[, 6], x$delsgs[, 5]),
           x = rep(1, nrow(x$delsgs))
  ) 
  return(res)
}

.obtain.weight <- function(
     coords = NULL,
     weight = NULL,
     weight.method = c("voronoi", "knn", "none"),
     ...
  ){
  params <- list(...)


  if (length(weight.method) > 1){
      weight.method <- match.arg(weight.method)
  }

  if (is.null(weight) && (weight.method %in% c("none", "voronoi", "knn"))){
      #if (is.integer(coords)) coords <- coords * 1.0
      #weight.mat <- pairDist(coords, coords)
      if (weight.method == 'knn'){
          if ("k" %in% names(params) && is.numeric(params$k)){
              k <- round(params$k)
          }else{
              k <- 10
          }
          #weight.mat <- .build.adj.m(weight.mat, k + 1) |> Matrix::t()
          #Matrix::diag(weight.mat) <- 0
          weight.mat <- .build.knn.adj(coords, k, weighted.distance = FALSE, bycol=FALSE)
      }else if (weight.method == "voronoi"){
          rlang::check_installed("deldir", "is required when `weight.method=='voronoi'`.") 
          weight.mat <- do.call(".build.adj.using.voronoi", list(coords, params)) 
      }
      weight.mat <- weight.mat * (1/Matrix::rowSums(weight.mat))
  }
  
  if (!is.null(weight)){
      if (inherits(weight, "listw") || inherits(weight, "nb")){
          weight.mat <- .convert_to_distmt(weight)
      }else if (inherits(weight, "matrix") && identical(rownames(weight), colnames(weight))){
          weight.mat <- weight |> Matrix::Matrix(sparse=TRUE)
      }else{
          cli::cli_abort(
            c("The {.var weight} should be a list (with named by sample_id) object containing `listw`, `nb`",
              "or squared `matrix` with equal the dimnames (the same to colnames of `data`.), Or it can be ",
              "a single `listw`, `nb` or `matrix` object, if it is provided (not NULL).")
          )
      }
  }

  if (is.null(weight) && !weight.method %in% c("none", "voronoi", "knn")){
      rlang::check_installed("spdep", paste0("is required to identify SVG or SVP with `", weight.method,"` ."))
      weight.mat <- do.call(weight.method, c(list(coords), params))
      weight.mat <- .convert_to_distmt(weight.mat)
  }

  return(weight.mat)
}

#' @importFrom deldir deldir
.build.adj.using.voronoi <- function(coords, flag.z = FALSE, ...){
  x <- coords[, 1]
  y <- coords[, 2]
  z <- NULL
  if (ncol(coords)>=3 && flag.z){
      z <- coords[,3]
  }
  
  res <- deldir::deldir(x, y, z, ...)
  res <- .convert_to_distmt(res)
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
                          random.seed = 123,
                          ...
                          ){

  rlang::check_installed(c("withr", "ks"), "is required to reproducible in the identification of SVG or SVP.")

  coords <- .normalize.coords(coords)
  
  lims <- c(range(coords[,1]), range(coords[,2]))
  h <- c(ks::hpi(coords[,1]), ks::hpi(coords[,2]))
 
  res <- withr::with_seed(random.seed, CalSpatialKldCpp(coords, x, lims, h, n, permutation))

  rownames(res) <- rownames(x)
  colnames(res) <- c("sp.kld", "boot.sp.kld.mean", "boot.sp.kld.sd", "pvalue")
  res <- cbind(res,
               padj = p.adjust(res[, "pvalue"], method = p.adjust.method)
            ) |> data.frame()
  res <- res |> dplyr::arrange(.data$padj, dplyr::desc(.data$sp.kld)) |>
         dplyr::mutate(rank = seq(nrow(res)))
  res <- res[match(rownames(x), rownames(res)), ]
  return(res)
}

.normalize.coords <- function(x){
  range_all <- max(apply(x, 2, function(col) diff(range(col))))
  x <- apply(x, 2, function(col) (col - min(col)) / range_all)
}

.weighted_by_hgt <- function(x, y){
  keep.names <- intersect(rownames(x), rownames(y))
  x <- MatElemMultiMat(x[keep.names,,drop=FALSE], y[keep.names,,drop=FALSE])
  rownames(x) <- keep.names
  colnames(x) <- colnames(y)
  return(x)
}
