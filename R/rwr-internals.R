#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom Rcpp evalCpp
.run_rwr <- function(g, 
                     edge.attr = 'weight',
                     seeds = NULL,
                     normalize.adj.method = c("laplacian", "row", "column", "none"),
                     restart = .7,
                     threads = 2L,
                     normalize.affinity = FALSE,
                     verbose = TRUE,
                     ...){

  normalize.adj.method <- match.arg(normalize.adj.method)
  if (inherits(g, 'igraph')){
    adj.m <- .extract.adj.m(g, edge.attr, verbose)
  }else{
    adj.m <- g
  }
  n.adj.m <- .normalize.adj.m(adj.m, normalize.adj.method)
  start.m <- .obtain.seeds(g, seeds = seeds, nrow = nrow(n.adj.m))
  
  if (restart > 1 || restart < 0){
    cli::cli_warn(c('The {.var restart} must be between 0 and 1.',
                    'it was set to .7 automatically.'))
    restart <- .7
  }

  stop.delta <- 1e-6
  stop.step <- 100
  tic()
  cli::cli_inform("Calculating the affinity score using random walk with restart ...")
  if (restart == 1){
    pt.m <- start.m
  }else{
    flag.threads <- Sys.getenv('RCPP_PARALLEL_NUM_THREADS')
    if (nchar(flag.threads)==0 && !is.null(threads)){
      RcppParallel::setThreadOptions(numThreads = threads)
    }
    pt.m <- parallelCalRWR(
                x = n.adj.m, 
                v = start.m, 
                restart = restart,
                stop_delta = stop.delta,
                stop_step = stop.step
            )
  }
  toc()

  pt.m[is.na(pt.m)] <- 0
  pt.m[pt.m < stop.delta] <- 0

  tic()
  cli::cli_inform(c('Tidying the result of random walk with restart ...'))
  if (ncol(pt.m) == 1){
    normalize.affinity <- FALSE
  }
  if (normalize.affinity){
    rlang::check_installed('broman', 'for `.run_rwr()` with `normalize.affinity = TRUE`.')
    pt.m <- broman::normalize(pt.m)
    pt.m[is.na(pt.m)] <- 0
  }
  if (ncol(pt.m) > 1){
    pt.m <- prop.table(pt.m, 1)
  }
  colnames(pt.m) <- colnames(start.m)
  rownames(pt.m) <- rownames(start.m)
  pt.m <- Matrix::Matrix(t(pt.m), sparse = TRUE)
  toc()
  return(pt.m)
}

.extract.adj.m <- function(g, edge.attr, verbose = FALSE){
  flag <- .check.edge.attr(g, edge.attr)
  if (inherits(flag, 'numeric')){
    adj.m <- igraph::as_adjacency_matrix(g, attr = edge.attr)
  }else{
    if (verbose){
      if (is.null(flag)){
        cli::cli_inform('The {.var edge.attr} does not exit in the {.var g}.')
      }else{
        cli::cli_inform('The {.var edge.attr} in the {.var g} is {.cls {class(flag)}}.')
      }
      cli::cli_inform(c('Unweighted adjacency matrix will be used.'))
    }
    adj.m <- igraph::as_adjacency_matrix(g)
  }
  return(adj.m)
}

.normalize.adj.m <- function(adj.m, normalize.adj.method){
  n.adj.m <- switch(
    normalize.adj.method,
    laplacian = {
      d <- Matrix::Diagonal(x=Matrix::colSums(adj.m)^-.5)
      d %*% adj.m %*% d
    },
    row = {
      d <- Matrix::Diagonal(x = Matrix::rowSums(adj.m)^-1)
      adj.m %*% d
    },
    column = {
      d <- Matrix::Diagonal(x = Matrix::colSums(adj.m)^-1)
      d %*% adj.m
    },
    none = adj.m,
  )
  return(n.adj.m)
}

.check.edge.attr <- function(g, edge.attr){
  x <- igraph::edge_attr(g, name='weight')
  return(x)
}

.obtain.seeds <- function(g, seeds = NULL, nrow){
  if (inherits(g, 'igraph')){    
    nm <- igraph::vertex_attr(g, name="name")
  }else{
    nm <- colnames(g)
  }
  if (is.null(seeds)){
    x <- Matrix::Matrix(diag(length(nm)), sparse = TRUE)
    rownames(x) <- colnames(x) <- nm
  }else{
    if (length(dim(seeds))!=2){
       cli::cli_abort(c('If {.var seeds} was provided, it must be a matrix',
                        'now it is a {.cls {class(seeds)}}'))
    }
    ind <- !is.na(match(rownames(seeds), nm))
    x <- Matrix::Matrix(0, nrow=nrow, ncol = ncol(seeds), sparse = TRUE)
    #x <- matrix(0, nrow = nrow, ncol = ncol(seeds))
    x[ind,] <- seeds[ind,]
    x <- x %*% diag(1/Matrix::colSums(x))
    x[is.na(x)] <- 0
    rownames(x) <- nm
    colnames(x) <- colnames(seeds)
    x <- x[,Matrix::colSums(x)>0,drop=FALSE]
    x <- Matrix::Matrix(x, sparse = TRUE)
  }
  return(x)
}

#calRWR <- function(x,
#                     v,
#                     restart = .75,
#                    delta = 1,
#                    step = 0,
#                     stop_delta = 1e-6, 
#                     stop_step = 50){
#    pt <- v
#    while(delta > stop_delta && step <= stop_step){
#      px <- (1 - restart) * x %*% pt + restart * v 
#      delta <- sum(abs(px - pt))
#      pt <- px
#      step <- step + 1
#    }
#    return(as.matrix(pt))
#}

