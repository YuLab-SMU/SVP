#' Calculate the activity of gene sets in spatial or single-cell data with restart walk with restart.
#' @description 
#' First, we calculated the distance between cells and between genes, between cells and genes in space
#' of \code{MCA}. Because the closer gene is to a cell, the more specific to such the cell it can 
#' be considered in \code{MCA} space (first reference). We extract the top nearest genes for each cells, 
#' to obtain the cells and cells association, genes and gens association, we also extract the top nearest 
#' cells or genes respectively, then combine all the association into the same network to obtain the adjacency 
#' matrix of all cells and genes. Another method is that we build the network using the combined \code{MCA} 
#' space of cells and genes directly, but this method can not combine the spatial physical space.(see also 
#' details). Next, we build a starting seed matrix (which each column measures the initial probability 
#' distribution of each gene set in graph nodes) for random walk with restart using the gene set and all
#' nodes of the graph, Finally we use restart walk with restart to calculate the affinity score for each
#' gene set or pathway.
#'
#' @rdname runSGSA-method
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
#' @param knn.consider.spcoord logical whether consider the spatial coordinates to run MCA, default is FALSE.
#' @param sp.alpha.add.weight only work when \code{knn.consider.spcoord=TRUE} and \code{knn.combined.cell.feature=FALSE},
#' which is weight of spatial space of the additive term in single cell and spatial space funsion formula, default is 0.2.
#' @param sp.beta.add.mp.weight only work when \code{knn.consider.spcoord=TRUE} and \code{knn.combined.cell.feature=FALSE},
#' which is weight of spatial space of the additive term and multiplicative term in single cell and spatial space funsion
#' formula, default is 0.1.
#' @param knn.used.reduction.dims the top components of the reduction with \code{MCA} to be used to build KNN 
#' Graph, default is 30.
#' @param knn.combined.cell.feature whether combined the embeddings of cells and features to find the nearest 
#' neighbor and build graph, default is FALSE, meaning the nearest neighbor will be found in cells to cells, 
#' features to features, cells to features respectively to build graph.
#' @param knn.graph.weighted logical whether consider the distance of nodes in the Nearest Neighbors, default is TRUE.
#' @param knn.k.use numeric the number of the Nearest Neighbors nodes, default is 600.
#' @param rwr.restart  default is 0.75.
#' @param rwr.normalize.adj.method character the method to normalize the adjacency matrix of the input graph,
#' default is \code{laplacian}.
#' @param rwr.normalize.affinity logical whether normalize the activity (affinity) result score using quantile normalisation,
#' default is FALSE.
#' @param rwr.prop.normalize logical whether divide the specific activity score by total activity score for a sample,
#' default is FALSE. 
#' @param rwr.threads the threads to run Random Walk With Restart (RWR), default is 2L.
#' @param hyper.test.weighted character which method to weight the activity score of cell, should is one of "Hypergeometric", "Wallenius", 
#' "none", default is "Hypergeometric".
#' @param hyper.test.by.expr logical whether using the expression matrix to find the nearest genes of cells, default is \code{TRUE},
#' if it is \code{FALSE}, meaning using the result of reduction to find the nearest genes of cells to perfrom the \code{hyper.test.weighted}.
#' @param cells Vector specifying the subset of cells to be used for the calculation of the activaty score or identification 
#' of SV features. This can be a character vector of cell names, an integer vector of column indices or a logical vector, 
#' default is NULL, meaning all cells to be used for the calculation of the activaty score or identification of SV features. 
#' @param features Vector specifying the subset of features to be used for the calculation of the activaty score or identification
#' of SV features. This can be a character vector of features names, an integer vector of row indices or a logical vector,
#' default is NULL, meaning all features to be used for the calculation of the activaty score or identification of SV features.
#' @param verbose logical whether print the intermediate message when running the program, default is TRUE.
#' @param ... additional parameters
#' @return a \linkS4class{SVPExperiment} or a \linkS4class{SingleCellExperiment}, see details.
#'
#' @details
#' if input is a \linkS4class{SVPExperiment}, output will be also a \linkS4class{SVPExperiment}, the activity score of gene sets
#' was stored in \code{assay} slot of the specified \code{gsvaexp}, and the spatially variable gene sets result is stored in \code{svDfs}
#' of the specified \code{gsvaexp}, which is a \linkS4class{SingleCellExperiment}. If input is a \linkS4class{SingleCellExperiment}
#' (which is extracted from \linkS4class{SVPExperiment} using \code{gsvaExp()} funtion), output will be also a
#' \linkS4class{SingleCellExperiment}, the activity score of gene sets result can be extracted using \code{assay()} function. The 
#' spatially variable gene sets result can be extracted using \code{svDf()} function.
#'
#' When the \code{knn.consider.spcoord = TRUE}, \code{combined.cell.feature=FALSE} and the input \code{data} contains the spatial space.
#' The distance between cells will be reconstructed by taking into account both the space of \code{MCA} from cell transcriptomics data and 
#' the physical space of cells in the following way (refer to the second refercence article):
#'
#' \eqn{C.dist = (1 - \beta) * ((1-\alpha) * S.dist + \alpha * P.dist) + \beta * S.dist \odot P.dist}
#' 
#' where \eqn{C.dist} is the new distance matrix of cells, \eqn{S.dist} is the distance matrix of cells in the \code{MCA} space from 
#' transcriptomics data, \eqn{P.dist} is the distance matrix from physical space of cells, \eqn{beta} weights the contributions of the 
#' additive and multiplicative terms, which is the argument \code{sp.beta.add.mp.weight}, \eqn{alpha} weighs the contributions of \eqn{P.dist} 
#' and \eqn{S.dist}, which is the argument \code{sp.alpha.add.weight}, and the \eqn{\odot} is the element-wise product.
#'
#' The affinity score is calculated in the following way (refer to the third refercence article):
#' 
#' \eqn{P_{t+1} = (1 - r) * M * P_{t} + r * P_{0}}
#' 
#' where \eqn{P_{0}} is the initial probability distribution for each gene set, \eqn{M} is the transition matrix that is the column normalization 
#' of adjacency matrix of graph, \eqn{r} is the the global restart probability, \eqn{P_{t+1}} and \eqn{P_{t}} are the probability distribution in 
#' each iteration. After several iterations, the difference between \eqn{P_{t+1}} and \eqn{P_{t}} becomes negligible, the stationary probability 
#' distribution is reached, and the elements for each gene set represent a proximity measure from every graph node. Iterations are stoped when the 
#' difference between \eqn{P_{t+1}} and \eqn{P_{t}} falls below 1e-6.
#' 
#' @references
#' 1. Cortal, A., Martignetti, L., Six, E. et al. Gene signature extraction and cell identity recognition at the single-cell 
#'    level with Cell-ID. Nat Biotechnol 39, 1095–1102 (2021). https://doi.org/10.1038/s41587-021-00896-6
#'
#' 2. Arutyunyan, A., Roberts, K., Troulé, K. et al. Spatial multiomics map of trophoblast development in early pregnancy. 
#'    Nature, 616, 143–151 (2023). https://doi.org/10.1038/s41586-023-05869-0.
#'
#' 3. Alberto Valdeolivas, Laurent Tichit, Claire Navarro, Sophie Perrin, et al. Random walk with restart on multiplex and 
#'    heterogeneous biological networks, Bioinformatics, 35, 3, 497–505(2019), https://doi.org/10.1093/bioinformatics/bty637
#'
#' @seealso [`cluster.assign`] to classify cell using the activity score of gene sets base \code{kmean} and [`runKldSVG`] to identify the 
#' spatiall variable or specified cell gene sets or a features.
#' @author Shuangbin Xu
#' @export
setGeneric('runSGSA', 
  function(
    data, 
    gset.idx.list,
    gsvaExp.name = 'gset1.rwr',
    min.sz = 10,
    max.sz = Inf,
    gene.occurrence.rate = .4,
    assay.type = 'logcounts',
    knn.consider.spcoord = FALSE,
    sp.alpha.add.weight = .2,
    sp.beta.add.mp.weight = .1,
    knn.used.reduction.dims = 30,
    knn.combined.cell.feature = FALSE,
    knn.graph.weighted = TRUE,
    knn.k.use = 600,
    rwr.restart = .75,
    rwr.normalize.adj.method = c("laplacian", "row", "column", "none"),
    rwr.normalize.affinity = FALSE,
    rwr.prop.normalize = FALSE,
    rwr.threads = 2L,
    hyper.test.weighted = c("Hypergeometric", "Wallenius", "none"),
    hyper.test.by.expr = TRUE,
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
    min.sz = 10, 
    max.sz = Inf,
    gene.occurrence.rate = .4,
    assay.type = 'logcounts',
    knn.consider.spcoord = FALSE,
    sp.alpha.add.weight = .2,
    sp.beta.add.mp.weight = .1,    
    knn.used.reduction.dims = 30,
    knn.combined.cell.feature = FALSE,
    knn.graph.weighted = TRUE,
    knn.k.use = 600,
    rwr.restart = .75,
    rwr.normalize.adj.method = c("laplacian", "row", "column", "none"),
    rwr.normalize.affinity = FALSE,
    rwr.prop.normalize = FALSE,
    rwr.threads = 2L,
    hyper.test.weighted = c("Hypergeometric", "Wallenius", "none"),
    hyper.test.by.expr = TRUE,
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
                     #consider.spcoord = knn.consider.spcoord, 
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

  gset.num <- .filter.gset.gene(features, gset.idx.list, min.sz, max.sz, gene.occurrence.rate)
  gset.idx.list <- gset.idx.list[match(rownames(gset.num), names(gset.idx.list))]

  flag1 <- .check_element_obj(data, key='spatialCoords', basefun=int_colData, namefun = names)
  if (flag1){
    coords <- .extract_element_object(data, key = 'spatialCoords', basefun=int_colData, namefun = names)
    coords <- .normalize.coords(coords)
  }else{
    coords <- NULL
  }  
  
  tic()
  cli::cli_inform(c("Building the nearest neighbor graph with the distance between 
                    features and cells ..."))
  
  rd.knn.gh <- .build.nndist.graph(
                       cells.rd,
                       features.rd,
                       coords,
                       knn.consider.spcoord,
                       sp.alpha.add.weight,
                       sp.beta.add.mp.weight,
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
  
  #gset.score.cells <- gset.score.cells[Matrix::rowSums(gset.score.cells) > 0,]
  features.expr <- assay(data, assay.type)
  
  if (hyper.test.weighted != 'none'){ 
      #cells.gene.num <- colSums(as.matrix(features.expr) > 0)
      if (hyper.test.by.expr){
          knn.gh <- .build.adj.m_by_expr(features.expr, top.n = knn.k.use, weighted.distance = knn.graph.weighted)
      }else{
          knn.gh <- rd.knn.gh[features, cells]
      }
      gset.hgt <- suppressWarnings(.run_hgt(
                           knn.gh,
                           seedstart.m[features,],
                           gset.idx.list,
                           #cells.gene.num,
                           m = gset.num[,'exp.gene.num'],
                           top.n = knn.k.use,
                           combined.cell.feature = knn.combined.cell.feature,
                           weighted.distance = knn.graph.weighted,
                           method = hyper.test.weighted
                  ))
      gset.score.cells <- gset.score.cells * gset.hgt
      #if (rwr.prop.normalize && ncol(gset.score.cells) > 1){
      #    gset.score.cells <- gset.score.cells |> prop.table(1) 
      #    gset.score.cells[is.na(gset.score.cells)] <- 0
      #}
  }
  
  gset.score.cells <- Matrix::Matrix(gset.score.cells, sparse = TRUE)
  gset.score.features <- .extract.features.rank(
                             gset.score.cells,
                             features.expr,
                             features,
                             gset.idx.list
                          )  
 
  x <- SingleCellExperiment(assays = list(affi.score = gset.score.cells))
  rowData(x) <- gset.num[rownames(x), ,drop=FALSE]
  
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