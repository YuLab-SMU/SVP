#' Calculate the activity of gene sets in spatial or single-cell data with restart walk with restart
#' and hyper test weighted.
#' @description 
#' First, we calculated the distance between cells and between genes, between cells and genes in space
#' of \code{MCA}. Because the closer gene is to a cell, the more specific to such the cell it can 
#' be considered in \code{MCA} space (first reference). We extract the top nearest genes for each cells, 
#' to obtain the cells and cells association, genes and gens association, we also extract the top nearest 
#' cells or genes respectively, then combine all the association into the same network to obtain the adjacency 
#' matrix of all cells and genes. Another method is that we build the network using the combined \code{MCA} 
#' space of cells and genes directly. Next, we build a starting seed matrix (which each column measures 
#' the initial probability distribution of each gene set in graph nodes) for random walk with restart using 
#' the gene set and all nodes of the graph. Finally, we employ the restart walk with restart algorithm to 
#' compute the affinity score for each gene set or pathway, which is then further weighted using the hypergeometric 
#' test result from the original expression matrix controlled by \code{hyper.test.weighted} parameter.
#'
#' @rdname runSGSA-method
#' @param data a \linkS4class{SingleCellExperiment} object normalized and have the result of 
#' \code{UMAP} or \code{TSNE}. Or a \linkS4class{SVPExperiment} object.
#' @param gset.idx.list gene set list contains the names, or GSON object or a gmt file, and
#' the online gmt file is also supported.
#' @param gsvaExp.name a character the name of \code{gsvaExp} of result \code{SVP} object.
#' @param symbol.from.gson logical whether extract the SYMBOL ID as \code{gset.idx.list}, only work
#' when \code{gset.idx.list} is a \code{GSON} object.
#' @param min.sz integer the minimum gene set number, default is 5, the number of gene sets 
#' smaller than \code{min.sz} will be ignored.
#' @param max.sz integer the maximum gene set number, default is Inf, the number of gene sets
#' larger than \code{max.sz} will be ignored.
#' @param gene.occurrence.rate the occurrence proportion of the gene set in the input object,
#' default is 0.2.
#' @param assay.type which expressed data to be pulled to build KNN Graph, default is \code{logcounts}.
#' @param knn.used.reduction.dims the top components of the reduction with \code{MCA} to be used to build KNN 
#' Graph, default is 30.
#' @param knn.combined.cell.feature whether combined the embeddings of cells and features to find the nearest 
#' neighbor and build graph, default is FALSE, meaning the nearest neighbor will be found in cells to cells, 
#' features to features, cells to features respectively to build graph.
#' @param knn.graph.weighted logical whether consider the distance of nodes in the Nearest Neighbors, default is TRUE.
#' @param knn.k.use numeric the number of the Nearest Neighbors nodes, default is 600.
#' @param rwr.restart the restart probability used for restart walk with restart, should be between 0 and 1, default is 0.75.
#' @param rwr.normalize.adj.method character the method to normalize the adjacency matrix of the input graph,
#' default is \code{laplacian}.
#' @param rwr.normalize.affinity logical whether normalize the activity (affinity) result score using quantile normalization,
#' default is FALSE.
#' @param rwr.prop.normalize logical whether divide the specific activity score by total activity score for a sample,
#' default is FALSE. 
#' @param rwr.threads the threads to run Random Walk With Restart (RWR), default is NULL, which will initialize with the default 
#' number of threads, you can also set this using \code{RcppParallel::setThreadOptions(numThreads=10)}.
#' @param hyper.test.weighted character which method to weight the activity score of cell, should is one of "Hypergeometric", "Wallenius", 
#' "none", default is "Hypergeometric".
#' @param hyper.test.by.expr logical whether using the expression matrix to find the nearest genes of cells, default is TRUE,
#' if it is FALSE, meaning using the result of reduction to find the nearest genes of cells to perform the \code{hyper.test.weighted}.
#' @param prop.score logical whether to normalize each feature for each sample, default is FALSE.
#' @param add.weighted.metric logical whether return the weight activity score of cell using the corresponding \code{hyper.test.weighted},
#' default is FALSE.
#' @param add.cor.features logical whether calculate the correlation between the new features and original features (genes), default
#' is FALSE. If it is TRUE the correlation result will be kept in fscoreDf which can be extracted using \code{fscoreDf()} function.
#' @param cells Vector specifying the subset of cells to be used for the calculation of the active score or identification 
#' of SV features. This can be a character vector of cell names, an integer vector of column indices or a logical vector, 
#' default is NULL, meaning all cells to be used for the calculation of the activity score or identification of SV features. 
#' @param features Vector specifying the subset of features to be used for the calculation of the activity score or identification
#' of SV features. This can be a character vector of features names, an integer vector of row indices or a logical vector,
#' default is NULL, meaning all features to be used for the calculation of the activity score or identification of SV features.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param ... additional parameters
#' @return a \linkS4class{SVPExperiment} or a \linkS4class{SingleCellExperiment}, see details.
#'
#' @details
#' if input is a \linkS4class{SVPExperiment}, output will be also a \linkS4class{SVPExperiment}, the activity score of gene sets
#' was stored in \code{assay} slot of the specified \code{gsvaexp}, and the spatially variable gene sets result is stored in \code{svDfs}
#' of the specified \code{gsvaexp}, which is a \linkS4class{SingleCellExperiment}. If input is a \linkS4class{SingleCellExperiment}
#' (which is extracted from \linkS4class{SVPExperiment} using \code{gsvaExp()} function), output will be also a
#' \linkS4class{SingleCellExperiment}, the activity score of gene sets result can be extracted using \code{assay()} function. The 
#' spatially variable gene sets result can be extracted using \code{svDf()} function.
#' The affinity score is calculated in the following way (refer to the second article):
#' 
#' \eqn{P_{t+1} = (1 - r) * M * P_{t} + r * P_{0}}
#' 
#' where \eqn{P_{0}} is the initial probability distribution for each gene set, \eqn{M} is the transition matrix that is the column normalization 
#' of adjacency matrix of graph, \eqn{r} is the global restart probability, \eqn{P_{t+1}} and \eqn{P_{t}} represent the probability distribution in 
#' each iteration. After several iterations, the difference between \eqn{P_{t+1}} and \eqn{P_{t}} becomes negligible, the stationary probability 
#' distribution is reached, indicating proximity measures from every graph node. Iterations are stopped when the difference between \eqn{P_{t+1}} 
#' and \eqn{P_{t}} falls below 1e-10.
#' 
#' @references
#' 1. Cortal, A., Martignetti, L., Six, E. et al. Gene signature extraction and cell identity recognition at the single-cell 
#'    level with Cell-ID. Nat Biotechnol 39, 1095–1102 (2021). https://doi.org/10.1038/s41587-021-00896-6
#'
#' 2. Alberto Valdeolivas, Laurent Tichit, Claire Navarro, Sophie Perrin, et al. Random walk with restart on multiplex and 
#'    heterogeneous biological networks, Bioinformatics, 35, 3, 497–505(2019), https://doi.org/10.1093/bioinformatics/bty637
#'
#' @seealso [`runDetectSVG`] and [`runKldSVG`] to identify the spatial variable features.
#' [`runGLOBALBV`] to explore the spatial co-distribution between the spatial variable features
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(sceSubPbmc)
#' library(SingleCellExperiment) |> suppressPackageStartupMessages()
#' library(scuttle) |> suppressPackageStartupMessages()
#' sceSubPbmc <- scuttle::logNormCounts(sceSubPbmc)
#' # the using runMCA to perform MCA (Multiple Correspondence Analysis)
#' # this is refer to the CelliD, but we using the Eigen to speed up.
#' # You can view the help information of runMCA using ?runMCA.
#' sceSubPbmc <- runMCA(sceSubPbmc, assay.type = 'logcounts')
#'
#' # Next, we can calculate the activity score of gene sets provided.
#' # Here, we use the Cell Cycle gene set from the Seurat 
#' # You can use other gene set, such as KEGG pathway, GO, Hallmark of MSigDB
#' # or TFs gene sets etc.
#' #
#' # supporting the list with names or gson object or the gmt file
#' # online gmt file is also be supported 
#' # such as 
#' # https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt
#' 
#' data(CellCycle.Hs)
#' sceSubPbmc <- runSGSA(sceSubPbmc, gset.idx.list = CellCycle.Hs, gsvaExp.name = 'CellCycle')
#' # Then a SVPE class which inherits SingleCellExperiment, is return.
#' sceSubPbmc
#' 
#' # You can obtaion the score matrix by following the commond
#' sceSubPbmc |> gsvaExp('CellCycle') 
#' sceSubPbmc |> gsvaExp("CellCycle") |> assay() |> t() |> head()
#' 
#' # Then you can use the ggsc or other package to visulize
#' # and you can try to use the findMarkers of scran or other packages to identify
#' # the different gene sets.   
#' #\donttest{
#'   library(ggplot2)
#'   library(ggsc)
#'   sceSubPbmc <- sceSubPbmc |> 
#'                 scater::runPCA(assay.type = 'logcounts', ntop = 600) |>
#'                 scater::runUMAP(dimred = 'PCA')
#'   # withReducedDim = TRUE, the original reducetion results from original gene features
#'   # will be add the colData in the sce.cellcycle.
#'   sce.cellcycle <- sceSubPbmc |> gsvaExp('CellCycle', withReducedDim=TRUE)
#'   sce.cellcycle
#'   sce.cellcycle |> sc_violin(
#'                       features = rownames(sce.cellcycle), 
#'                       mapping = aes(x=seurat_annotations, fill = seurat_annotations)
#'                    ) + 
#'                    scale_x_discrete(guide=guide_axis(angle=-45))
#'   sce.cellcycle |> sc_feature(features= "S", reduction='UMAP')
#'   library(scran)
#'   cellcycle.test.res <- sce.cellcycle |> findMarkers(
#'                      group = sce.cellcycle$seurat_annotations, 
#'                      test.type = 'wilcox', 
#'                      assay.type = 'affi.score', 
#'                      add.summary = TRUE
#'                   )
#'   cellcycle.test.res$B
#' #}
setGeneric('runSGSA', 
  function(
    data, 
    gset.idx.list,
    gsvaExp.name = 'gset1.rwr',
    symbol.from.gson = FALSE,
    min.sz = 5,
    max.sz = Inf,
    gene.occurrence.rate = .2,
    assay.type = 'logcounts',
    knn.used.reduction.dims = 30,
    knn.combined.cell.feature = FALSE,
    knn.graph.weighted = TRUE,
    knn.k.use = 600,
    rwr.restart = .75,
    rwr.normalize.adj.method = c("laplacian", "row", "column", "none"),
    rwr.normalize.affinity = FALSE,
    rwr.prop.normalize = FALSE,
    rwr.threads = NULL,
    hyper.test.weighted = c("Hypergeometric", "Wallenius", "none"),
    hyper.test.by.expr = TRUE,
    prop.score = FALSE,
    add.weighted.metric = FALSE,
    add.cor.features = FALSE,
    cells = NULL,
    features = NULL,
    verbose = TRUE, 
    ...
  )
  standardGeneric('runSGSA')
)

#' @importFrom SingleCellExperiment reducedDim<- reducedDimNames SingleCellExperiment
#' @importFrom SummarizedExperiment rowData colData<- rowData<- 
#' @importFrom pracma tic toc
#' @rdname runSGSA-method
#' @aliases runSGSA,SingleCellExperiment
#' @export runSGSA
setMethod('runSGSA', 
  'SingleCellExperiment',
  function(
    data,
    gset.idx.list,
    gsvaExp.name = 'gset1.rwr',
    symbol.from.gson = FALSE,
    min.sz = 5, 
    max.sz = Inf,
    gene.occurrence.rate = .2,
    assay.type = 'logcounts',
    knn.used.reduction.dims = 30,
    knn.combined.cell.feature = FALSE,
    knn.graph.weighted = TRUE,
    knn.k.use = 600,
    rwr.restart = .75,
    rwr.normalize.adj.method = c("laplacian", "row", "column", "none"),
    rwr.normalize.affinity = FALSE,
    rwr.prop.normalize = FALSE,
    rwr.threads = NULL,
    hyper.test.weighted = c("Hypergeometric", "Wallenius", "none"),
    hyper.test.by.expr = TRUE,
    prop.score = FALSE,
    add.weighted.metric = FALSE,
    add.cor.features = FALSE,
    cells = NULL,
    features = NULL,
    verbose = TRUE,
    ...
  ){
  hyper.test.weighted <- match.arg(hyper.test.weighted)
  knn.used.reduction <- 'MCA'
  if (!"MCA" %in% reducedDimNames(data)){
      cli::cli_warn(c("The {.cls {class(data)}} does not have MCA, run 'runMCA()' first."))
      data <- runMCA(data, 
                     assay.type = assay.type, 
                     ncomponents = knn.used.reduction.dims, 
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
  gset.idx.list <- .extract.gset(gset.idx.list, gene.name = symbol.from.gson)
  gset.num <- .filter.gset.gene(features, gset.idx.list, min.sz, max.sz, gene.occurrence.rate)
  gset.idx.list <- gset.idx.list[match(rownames(gset.num), names(gset.idx.list))]

  
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
  cli::cli_inform("Building the seed matrix using the gene set and the nearest neighbor 
                   graph for random walk with restart ...")

  seedstart.m <- .generate.gset.seed(rd.knn.gh, gset.idx.list)
  toc()

  gset.score <- .run_rwr(
                  rd.knn.gh,
                  cells, 
                  edge.attr = 'weight',
                  seeds = seedstart.m,
                  normalize.adj.method = rwr.normalize.adj.method,
                  restart = rwr.restart,
                  threads = rwr.threads,
                  normalize.affinity = rwr.normalize.affinity,
                  prop.normalize = rwr.prop.normalize,
                  verbose = verbose
                )

  gset.score.cells <- gset.score[, cells, drop=FALSE]
  
  features.expr <- assay(data, assay.type)
  features.expr <- features.expr[features, cells] 
  if (hyper.test.weighted != 'none'){ 
      if (hyper.test.by.expr){
          if (hyper.test.weighted == 'Hypergeometric') knn.graph.weighted <- FALSE
          knn.gh <- .build_col_knn(features.expr, knn.k.use, knn.graph.weighted)
      }else{
          knn.gh <- rd.knn.gh[features, cells]
      }
      gset.hgt <- suppressWarnings(.run_hgt(
                           fs2cell.adj = knn.gh,
                           fs2gset.adj = seedstart.m[features,],
                           gene.set.list = gset.idx.list[rownames(gset.score.cells)],
                           m = gset.num[,'exp.gene.num'],
                           top.n = knn.k.use,
                           combined.cell.feature = knn.combined.cell.feature,
                           weighted.distance = knn.graph.weighted,
                           method = hyper.test.weighted
                  ))
      gset.score.cells2 <- .weighted_by_hgt(gset.score.cells, gset.hgt)
      gset.score.cells2 <- .normalize_score(gset.score.cells2, prop.score)
      assay.res <- list(affi.score = as(gset.score.cells2, 'dgCMatrix'))
      if (add.weighted.metric){
          assay.res <- c(assay.res, 
                         list(rwr.score = as(gset.score.cells[rownames(gset.score.cells2),,drop=FALSE], 'dgCMatrix'),
                            hyper.weighted = as(gset.hgt,'dgCMatrix')))
      }
  }else{
      gset.score.cells <- .normalize_score(gset.score.cells, prop.score)
      assay.res <- list(affi.score = as(gset.score.cells, 'dgCMatrix'))
  }
  
  x <- SingleCellExperiment(assays = assay.res)
  rowData(x) <- gset.num[rownames(x), ,drop=FALSE]

  if (add.cor.features){
      gset.score.features <- .extract.features.rank(
                                 assay(x),
                                 features.expr,
                                 features,
                                 gset.idx.list
                              )
      
      x <- .add.int.rowdata(sce=x, getfun=fscoreDfs, 
                            setfun1 = `fscoreDfs<-`, 
                            setfun2 = `fscoreDf<-`, 
                            namestr = "rwr.score", 
                            val = gset.score.features)
  }
  da <- .sce_to_svpe(data) 
  gsvaExp(da, gsvaExp.name) <- x
  new.reduced <- .build.new.reduced(rd.df, cells, features, rd.f.nm)
  reducedDim(da, knn.used.reduction) <- new.reduced
  return(da)
})
