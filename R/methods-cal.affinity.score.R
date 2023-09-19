setGeneric('cal.affinity.score', 
  function(
    data, 
    gset.idx.list,
    reduction.name = 'MCA', 
    dims = 30, 
    knn.k.use = 350,
    restart = 0.75, 
    subset.row = NULL,
    ...
  )
  standardGeneric('cal.affinity.score')
)

setMethod('cal.affinity.score', 
  'SingleCellExperiment',
  function(
    data,
    gset.idx.list,
    reduction.name = 'MCA',
    knn.fun = 'findKmknn',
    knn.BPPARAM = SerialParam(),
    knn.graph.weighted = TRUE,
    knn.graph.directed = FALSE,
    dims = 30,
    knn.k.use = 350,
    restart = .7,
    cells = NULL,
    features = NULL,
    ...
  ){
  if (reduction.name %in% reducedDimNames(data)){
      cli::cli_abort(c("The {.cls {class(data)}} does not have ", reduction.name))
  }
  mca.res <- reducedDim(data, reduction.name)
  mca.res <- rbind(mca.res, attr(mca.res, 'featuresCoords'))

  mca.res <- .subset_data(x = mca.res, n = c(cells, features))
  
  dims <- min(ncol(mca.res), dims) 
  
  mca.res <- mca.res[, seq(dims), drop = FALSE]

  mca.knn.gh <- .build.knn.graph(
                  mca.res, 
                  knn.k.use = knn.k.use, 
                  fun.nm = knn.fun,
                  BPPARAM = knn.BPPARAM,
                  weighted.distance = knn.graph.weighted,
                  graph.directed = knn.graph.directed,
                  ...
                )

  gset.score <- .cal.gset.score(mca.knn.gh)

})
