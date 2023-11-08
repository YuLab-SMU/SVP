#' @title Calculate the activity of gene sets in spatial or single-cell data with Restart Walk with Restart 
#' @rdname sc.rwr-method
#' @param data a \linkS4class{SingleCellExperiment} object normalized and have the result of 
#' \code{UMAP} or \code{TSNE}. Or a \linkS4class{SVPExperiment} object.
#' @param gset.idx.list gene set list contains the names.
#' @param gsvaExp.name a character the name of \code{gsvaExp} of result \code{SVP} object.
#' @param min.sz integer the minimum gene set number, default is 10, the number of gene sets 
#' smaller than \code{min.sz} will be ignored.
#' @param max.sz integer the maximum gene set number, default is Inf, the number of gene sets
#' larger than \code{max.sz} will be ignored.
#' @param gene.occurrence.rate the occurrence proportion of the gene set in the input object,
#' default is 0.4.
#' @param assay.type which expressed data to be pulled to build KNN Graph, default is \code{logcounts}.
#' @param knn.mca.consider.spcoord logical whether consider the spatial coordinates to run MCA, default is TRUE.
#' @param knn.used.reduction.dims the top components of the reduction with \code{MCA} to be used to build KNN 
#' Graph, default is 30.
#' @param knn.combined.cell.feature whether combined the embeddings of cells and features to find the nearest 
#' neighbor and build graph, default is FALSE, meaning the nearest neighbor will be found in cells to cells, 
#' features to features, cells to features respectively to build graph.
#' @param knn.graph.weighted logical whether consider the distance of nodes in the Nearest Neighbors, default is TRUE.
#' @param knn.k.use numeric the number of the Nearest Neighbors nodes, default is 600.
#' @param rwr.restart  default is 0.7.
#' @param rwr.normalize.adj.method character the method to normalize the adjacency matrix of the input graph,
#' default is \code{laplacian}.
#' @param rwr.normalize.affinity logical whether normalize the activity (affinity) result score using quantile normalisation,
#' default is TRUE.
#' @param rwr.threads the threads to run Random Walk With Restart (RWR), default is 2L.
#' @param cells Vector specifying the subset of cells to be used for the calculation of the activaty score or identification 
#' of SV features. This can be a character vector of cell names, an integer vector of column indices or a logical vector, 
#' default is NULL, meaning all cells to be used for the calculation of the activaty score or identification of SV features. 
#' @param features Vector specifying the subset of features to be used for the calculation of the activaty score or identification
#' of SV features. This can be a character vector of features names, an integer vector of row indices or a logical vector,
#' default is NULL, meaning all features to be used for the calculation of the activaty score or identification of SV features.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param ... additional parameters
#' @export
setGeneric('sc.rwr', 
  function(
    data, 
    gset.idx.list,
    gsvaExp.name = 'gset1.rwr',
    min.sz = 10,
    max.sz = Inf,
    gene.occurrence.rate = .4,
    assay.type = 'logcounts',
    knn.mca.consider.spcoord = TRUE,
    knn.used.reduction.dims = 30,
    knn.combined.cell.feature = FALSE,
    knn.graph.weighted = TRUE,
    knn.k.use = 600,
    rwr.restart = .7,
    rwr.normalize.adj.method = c("laplacian","row","column","none"),
    rwr.normalize.affinity = TRUE,
    rwr.threads = 2L,
    cells = NULL,
    features = NULL,
    verbose = TRUE, 
    ...
  )
  standardGeneric('sc.rwr')
)

#' @importFrom SingleCellExperiment reducedDim<- reducedDimNames SingleCellExperiment
#' @importFrom SummarizedExperiment rowData colData<- rowData<- 
#' @importFrom pracma tic toc
#' @rdname sc.rwr-method
#' @aliases sc.rwr,SingleCellExperiment
#' @export sc.rwr
setMethod('sc.rwr', 
  'SingleCellExperiment',
  function(
    data,
    gset.idx.list,
    gsvaExp.name = 'gset1.rwr',
    min.sz = 10, 
    max.sz = Inf,
    gene.occurrence.rate = .4,
    assay.type = 'logcounts',
    knn.mca.consider.spcoord = TRUE,
    knn.used.reduction.dims = 30,
    knn.combined.cell.feature = FALSE,
    knn.graph.weighted = TRUE,
    knn.k.use = 600,
    rwr.restart = .7,
    rwr.normalize.adj.method = c("laplacian","row","column","none"),
    rwr.normalize.affinity = TRUE,
    rwr.threads = 2L,
    cells = NULL,
    features = NULL,
    verbose = TRUE,
    ...
  ){
  knn.used.reduction <- 'MCA'
  if (!"MCA" %in% reducedDimNames(data)){
      cli::cli_warn(c("The {.cls {class(data)}} does not have MCA, run 'runMCA()' first."))
      data <- runMCA(data, 
                     assay.type = assay.type, 
                     ncomponents = knn.used.reduction.dims, 
                     consider.spcoord = knn.mca.consider.spcoord, 
                     subset.row = cells, 
                     subset.col = features)
  }
  rd.df <- reducedDim(data, knn.used.reduction)
  rd.f.nm <- switch(knn.used.reduction, MCA='genesCoordinates', PCA='rotation')
  rd.f.res <- attr(rd.df, rd.f.nm)

  cells <- .subset_ind(rd.df, cells)
  features <- .subset_ind(rd.f.res, features)

  dims <- min(ncol(rd.df), knn.used.reduction.dims)
  cells.rd <- rd.df[cells, seq(dims), drop=FALSE]
  features.rd <- rd.f.res[features, seq(dims), drop=FALSE]

  gset.num <- .filter.gset.gene(features, gset.idx.list)
  gset.idx.list <- gset.idx.list[names(gset.idx.list) %in% rownames(gset.num)]

  tic()
  cli::cli_inform(c("Building the nearest neighbor graph with the distance between 
		    features and cells ..."))
  
  rd.knn.gh <- .build.nndist.graph(
                       cells.rd,
                       features.rd,
                       top.n = knn.k.use,
                       combined.cell.feature = knn.combined.cell.feature,
                       weighted.distance = knn.graph.weighted
               )
  toc()

  tic()
  cli::cli_inform("Building the seed matrix using the gene set and the neareast 
		   neighbor graph for Random Walk with Restart ...")

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

  da <- .sce_to_svpe(data) 
  gsvaExp(da, gsvaExp.name) <- x
  new.reduced <- .build.new.reduced(rd.df, cells, features, rd.f.nm)
  reducedDim(da, knn.used.reduction) <- new.reduced
  return(da)
})
