#' @export
setGeneric('cal.affinity.score', 
  function(
    data, 
    gset.idx.list,
    reduction.name = 'MCA', 
    dims = 30,
    knn.fun = findKmknn,
    knn.BPPARAM = SerialParam(),
    knn.graph.weighted = TRUE,
    knn.graph.directed = FALSE,
    knn.k.use = 500,
    rwr.restart = .7,
    rwr.normalize.adj.method = c("laplacian","row","column","none"),
    rwr.normalize.affinity = TRUE,
    rwr.BPPARAM = SerialParam(),
    rwr.verbose = TRUE,
    cells = NULL,
    features = NULL, 
    ...
  )
  standardGeneric('cal.affinity.score')
)

#' @importFrom pracma tic toc
#' @export
setMethod('cal.affinity.score', 
  'SingleCellExperiment',
  function(
    data,
    gset.idx.list,
    reduction.name = 'MCA',
    dims = 30,
    knn.fun = 'findKmknn',
    knn.BPPARAM = SerialParam(),
    knn.graph.weighted = TRUE,
    knn.graph.directed = FALSE,
    knn.k.use = 500,
    rwr.restart = .7,
    rwr.normalize.adj.method = c("laplacian","row","column","none"),
    rwr.normalize.affinity = TRUE,
    rwr.BPPARAM = SerialParam(),
    rwr.verbose = TRUE,
    cells = NULL,
    features = NULL,
    ...
  ){
  if (!reduction.name %in% reducedDimNames(data)){
      cli::cli_abort(c("The {.cls {class(data)}} does not have {.var reduction.name}"))
  }
  rd.res <- reducedDim(data, reduction.name)
  rd.f.res <- attr(rd.res, 'featuresCoords')

  cells <- .subset_ind(rd.res, cells)
  features <- .subset_ind(rd.f.res, features)

  rd.res <- rbind(rd.res, rd.f.res)
  rd.res <- .subset_data(x = rd.res, n = c(cells, features))
  
  dims <- min(ncol(rd.res), dims) 
  
  rd.res <- rd.res[, seq(dims), drop = FALSE]

  tic()
  cli::cli_inform(c("Building the knn graph ..."))

  rd.knn.gh <- .build.knn.graph(
                  rd.res, 
                  knn.k.use = knn.k.use, 
                  fun.nm = knn.fun,
                  BPPARAM = knn.BPPARAM,
                  weighted.distance = knn.graph.weighted,
                  graph.directed = knn.graph.directed,
                  ...
                )
  toc()

  tic()
  cli::cli_inform(c("Building the seed matrix using the gene set ",
		    "and the KNN graph built before for Random Walk ",
		    "with Restart ..."))

  seedstart.m <- .generate.gset.seed(rd.knn.gh, gset.idx.list)
  toc()

  tic()
  cli::cli_inform(c('Calculating the affinity score using Random Walk with Restart ...'))
  gset.score <- .run_rwr(
                  rd.knn.gh, 
		  edge.attr = 'weight',
                  seeds = seedstart.m,
                  normalize.adj.method = rwr.normalize.adj.method,
                  restart = rwr.restart,
		  BPPARAM = rwr.BPPARAM,
                  normalize.affinity = rwr.normalize.affinity,
                  verbose = rwr.verbose
                )
  gset.score <- gset.score[, cells]
  x <- SingleCellExperiment(assays = list(affi.score = gset.score))
  data <- .sce_to_svpe(data) 
  gsvaExp(data, "rwr") <- x
  toc()
  return(data)

})
