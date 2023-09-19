setGeneric('cal.affinity.score', 
  function(
    data, 
    gset.idx.list,
    reduction.name = 'MCA', 
    dims=20, 
    knn.k.use = 300,
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
    dims = 20,
    knn.k.use = 300,
    restart = .75,
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


})
