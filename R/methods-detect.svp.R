#' @title Calculate the activity of gene sets in spatial or single-cell data with Restart Walk 
#' with Restart and detect the spatially or single cell variable gene sets with Kullback–Leibler 
#' divergence of 2D Weighted Kernel Density
#' @rdname detect.svp-method
#' @param data a \linkS4class{SingleCellExperiment} object normalized and done \code{MCA} and 
#' \code{UMAP} or \code{TSNE}.
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
#' @param knn.used.reduction.dims the top components of the reduction with \code{knn.used.reduction} 
#' to be used to build KNN Graph, default is 30.
#' @param knn.BPPARAM A BiocParallelParam object specifying whether the built of KNN should be parallelized
#' default is \code{SerialParam()}, meaning no parallel. You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=T)}
#' to parallel it, the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}.
#' @param knn.graph.weighted logical whether consider the distance of nodes in the Nearest Neighbors, default is TRUE.
#' @param knn.graph.directed logical whether consider the direction of the KNN Graph, default is FALSE.
#' @param knn.k.use numeric the number of the Nearest Neighbors nodes, default is 500.
#' @param rwr.restart  default is 0.7.
#' @param rwr.normalize.adj.method character the method to normalize the adjacency matrix of the input graph,
#' default is \code{laplacian}.
#' @param rwr.normalize.affinity logical whether normalize the activity (affinity) result score using quantile normalisation,
#' default is TRUE.
#' @param rwr.threads the threads to run Random Walk With Restart (RWR), default is 2L.
#' @param sv.used.reduction character used as spatial coordinates to detect SVG, default is \code{UMAP},
#' if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates.
#' @param sv.grid.n numeric number of grid points in the two directions to estimate 2D weighted kernel density, default is 100.
#' @param sv.permutation numeric the number of permutation for each single feature to detect the signicantly spatially or single 
#' cell variable features, default is 100.
#' @param sv.p.adjust.method character the method to adjust the pvalue of the result, default is \code{bonferroni}.
#' @param sv.BPPARAM A BiocParallelParam object specifying whether the identification of SV features should be parallelized
#' default is \code{SerialParam()}, meaning no parallel. You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=T)}
#' to parallel it, the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}. default is \code{SerialParam()}.
#' @param run.sv logical whether run the identication of SV features using \code{kldSVG}, if it is \code{FALSE}, the identification 
#' of SV features will not be done, default is TRUE.
#' @param cells Vector specifying the subset of cells to be used for the calculation of the activaty score or identification 
#' of SV features. This can be a character vector of cell names, an integer vector of column indices or a logical vector, 
#' default is NULL, meaning all cells to be used for the calculation of the activaty score or identification of SV features. 
#' @param features Vector specifying the subset of features to be used for the calculation of the activaty score or identification
#' of SV features. This can be a character vector of features names, an integer vector of row indices or a logical vector,
#' default is NULL, meaning all features to be used for the calculation of the activaty score or identification of SV features.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param random.seed numeric random seed number to repeatability, default is 1024.
#' @param ... additional parameters
#' @export
setGeneric('detect.svp', 
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
    sv.p.adjust.method = "bonferroni",
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

#' @importFrom SingleCellExperiment reducedDim<- reducedDimNames SingleCellExperiment
#' @importFrom SummarizedExperiment rowData colData<- rowData<- 
#' @importFrom pracma tic toc
#' @importFrom BiocNeighbors findKmknn
#' @rdname detect.svp-method
#' @aliases detect.svp,SingleCellExperiment
#' @export detect.svp
setMethod('detect.svp', 
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
    sv.p.adjust.method = "bonferroni",
    sv.BPPARAM = SerialParam(),
    run.sv = TRUE,
    cells = NULL,
    features = NULL,
    verbose = TRUE,
    random.seed = 1024,
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

  #rd.res <- rbind(rd.df, rd.f.res)
  #rd.res <- .subset_data(x = rd.res, n = c(cells, features))
  #
  #dims <- min(ncol(rd.res), knn.used.reduction.dims) 
  #
  #rd.res <- rd.res[, seq(dims), drop = FALSE]
  
  dims <- min(ncol(rd.df), knn.used.reduction.dims)
  cells.rd <- rd.df[cells, seq(dims), drop=FALSE]
  features.rd <- rd.f.res[features, seq(dims), drop=FALSE]

  gset.num <- .filter.gset.gene(features, gset.idx.list)

  gset.idx.list <- gset.idx.list[names(gset.idx.list) %in% rownames(gset.num)]

  tic()
  cli::cli_inform(c("Building the knn graph ..."))

  #rd.knn.gh <- .build.knn.graph(
  #                rd.res, 
  #                knn.k.use = knn.k.use, 
  #                fun.nm = findKmknn,
  #                BPPARAM = knn.BPPARAM,
  #                weighted.distance = knn.graph.weighted,
  #                graph.directed = knn.graph.directed,
  #                ...
  #             )
  rd.knn.gh <- .build.nndist.graph(
                   cells.rd, 
                   features.rd, 
                   top.n = knn.k.use, 
                   weighted.distance = knn.graph.weighted, 
                   graph.directed = knn.graph.directed
               )
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
                      Kullback-Leibler divergence of 2D Weighted Kernel Density ...") 
      
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
