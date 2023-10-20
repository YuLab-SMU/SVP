#' @export
setGeneric('detect.svp', 
  function(
    data, 
    gset.idx.list,
    gsvaExp.name = 'gset1.rwr',
    dims = 30,
    min.sz = 10,
    max.sz = Inf,
    gene.occurrence.rate = .4,
    knn.used.reduction = c('MCA', 'PCA'),
    knn.used.reduction.dims = 30,    
    knn.BPPARAM = SerialParam(),
    knn.graph.weighted = TRUE,
    knn.graph.directed = FALSE,
    knn.k.use = 500,
    rwr.restart = .7,
    rwr.normalize.adj.method = c("laplacian","row","column","none"),
    rwr.normalize.affinity = TRUE,
    rwr.threads = 2L,
    sv.used.reduction = c('UMAP', 'TSNE'),
    sv.grid.n = 100,
    sv.permutation = 100,
    sv.p.adjust.method = "fdr",
    sv.BPPARAM = SerialParam(),
    run.sv = TRUE, 
    cells = NULL,
    features = NULL,
    verbose = TRUE, 
    random.seed = 1024,
    ...
  )
  standardGeneric('detect.svp')
)

#' @importFrom SingleCellExperiment rowData
#' @importFrom pracma tic toc
#' @importFrom BiocNeighbors findKmknn
#' @importFrom withr with_seed
#' @export
setMethod('detect.svp', 
  'SingleCellExperiment',
  function(
    data,
    gset.idx.list,
    gsvaExp.name = 'gset1.rwr',
    min.sz = 10, 
    max.sz = Inf,
    gene.occurrence.rate = .4,
    knn.used.reduction = c('MCA', 'PCA'),
    knn.used.reduction.dims = 30,
    knn.BPPARAM = SerialParam(),
    knn.graph.weighted = TRUE,
    knn.graph.directed = FALSE,
    knn.k.use = 500,
    rwr.restart = .7,
    rwr.normalize.adj.method = c("laplacian","row","column","none"),
    rwr.normalize.affinity = TRUE,
    rwr.threads = 2L,
    sv.used.reduction = c('UMAP', 'TSNE'),
    sv.grid.n = 100,
    sv.permutation = 100,
    sv.p.adjust.method = "fdr",
    sv.BPPARAM = SerialParam(),
    run.sv = TRUE,
    cells = NULL,
    features = NULL,
    verbose = TRUE,
    random.seed = 1024,
    ...
  ){

  knn.used.reduction <- match.arg(knn.used.reduction)
  runrd <- c('MCA'='runMCA()', 'PCA' = 'runPCA()')[knn.used.reduction]
  if (!knn.used.reduction %in% reducedDimNames(data)){
      cli::cli_abort(c("The {.cls {class(data)}} does not have knn.used.reduction = {knn.used.reduction},
		       please run { runrd} first."))
  }
  rd.df <- reducedDim(data, knn.used.reduction)
  rd.f.nm <- switch(knn.used.reduction, MCA='genesCoordinates', PCA='rotation')
  rd.f.res <- attr(rd.df, rd.f.nm)

  cells <- .subset_ind(rd.df, cells)
  features <- .subset_ind(rd.f.res, features)

  rd.res <- rbind(rd.df, rd.f.res)
  rd.res <- .subset_data(x = rd.res, n = c(cells, features))
  
  dims <- min(ncol(rd.res), dims) 
  
  rd.res <- rd.res[, seq(dims), drop = FALSE]

  gset.num <- .filter.gset.gene(features, gset.idx.list)

  gset.idx.list <- gset.idx.list[names(gset.idx.list) %in% rownames(gset.num)]

  tic()
  cli::cli_inform(c("Building the knn graph ..."))

  rd.knn.gh <- withr::with_seed(random.seed, .build.knn.graph(
                  rd.res, 
                  knn.k.use = knn.k.use, 
                  fun.nm = "findKmknn",
                  BPPARAM = knn.BPPARAM,
                  weighted.distance = knn.graph.weighted,
                  graph.directed = knn.graph.directed,
                  ...
               ))
  toc()

  tic()
  cli::cli_inform("Building the seed matrix using the gene set and the KNN graph 
		   built before for Random Walk with Restart ...")

  seedstart.m <- .generate.gset.seed(rd.knn.gh, gset.idx.list)
  toc()

  gset.score <- .run_rwr(
                  rd.knn.gh, 
                  edge.attr = 'weight',
                  seeds = seedstart.m,
                  normalize.adj.method = rwr.normalize.adj.method,
                  restart = rwr.restart,
                  threads = rwr.threads,
                  normalize.affinity = rwr.normalize.affinity,
                  verbose = verbose
                )

  gset.score.cells <- gset.score[, cells]

  gset.score.cells <- gset.score.cells[Matrix::rowSums(gset.score.cells) > 0,]

  gset.score.features <- .extract.features.score(
                            gset.score, 
                            rownames(gset.score.cells), 
                            features, 
                            gset.idx.list
                         )

  x <- SingleCellExperiment(assays = list(affi.score = gset.score.cells))

  rowData(x) <- gset.num[rownames(gset.num) %in% rownames(x), ]
  
  x <- .add.int.rowdata(sce=x, getfun=fscoreDfs, 
			setfun1 = `fscoreDfs<-`, 
			setfun2 = `fscoreDf<-`, 
			namestr = "rwr.score", 
			val = gset.score.features)

  flag1 <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)

  flag2 <- any(sv.used.reduction %in% reducedDimNames(data))

  if((flag1 || flag2) && run.sv){
      if (flag2){
          coords <- reducedDim(data, sv.used.reduction)
	  coords <- coords[,c(1, 2)]
      }
      if (flag1){
          coords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
      }
      tic()
      cli::cli_inform("Identifying the spatially variable gene sets (pathway) based on 
                      Kullback–Leibler divergence of 2D Weighted Kernel Density ...") 
      
      res.sv <- .identify.svg(
			gset.score.cells, 
                        coords = coords,
                        n = sv.grid.n,
                        permutation = sv.permutation,
			p.adjust.method = sv.p.adjust.method,
                        BPPARAM = sv.BPPARAM, 
                        random.seed = random.seed,
                        ...)
      
      x <- .add.int.rowdata(sce = x, 
			    getfun = svDfs, 
			    setfun1 = `svDfs<-`, 
			    setfun2 = `svDf<-`, 
			    namestr = 'sv.kld', 
                            val = res.sv)
      toc()
  }
  
  da <- .sce_to_svpe(data) 
  gsvaExp(da, gsvaExp.name) <- x
  new.reduced <- .build.new.reduced(rd.df, cells, features, rd.f.nm)
  reducedDim(da, knn.used.reduction) <- new.reduced
  return(da)
})
