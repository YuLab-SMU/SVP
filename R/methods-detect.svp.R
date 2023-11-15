#' @title detect.svp
#' @description Calculate the activity of gene sets in spatial or single-cell data with Restart Walk 
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
#' @param knn.consider.spcoord logical whether combine the space of \code{MCA} and the spatial physical 
#' space, default is FALSE. It only works when the input \code{data} has spatial coordinates and the
#' argument \code{knn.combined.cell.feature=FALSE} and \code{knn.consider.spcoord=TRUE}.
#' @param sp.alpha.add.weight only work when \code{knn.consider.spcoord=TRUE} and \code{knn.combined.cell.feature=FALSE},
#' which is weight of spatial space of the additive term in single cell and spatial space funsion formula, default is 0.2.
#' @param sp.beta.add.mp.weight only work when \code{knn.consider.spcoord=TRUE} and \code{knn.combined.cell.feature=FALSE},
#' which is weight of spatial space of the additive term and multiplicative term in single cell and spatial space funsion 
#' formula, default is 0.1.
#' @param knn.used.reduction.dims the top components of the reduction with \code{knn.used.reduction} 
#' to be used to build KNN Graph, default is 30.
#' @param knn.combined.cell.feature whether combined the embeddings of cells and features to find the nearest
#' neighbor and build graph, default is FALSE, meaning the nearest neighbor will be found in cells to cells,
#' features to features, cells to features respectively to build graph.
#' @param knn.graph.weighted logical whether consider the distance of nodes in the nearest neighbors, default is TRUE.
#' @param knn.k.use numeric the number of the Nearest Neighbors nodes, default is 400.
#' @param rwr.restart  default is 0.75.
#' @param rwr.normalize.adj.method character the method to normalize the adjacency matrix of the input graph,
#' default is \code{laplacian}.
#' @param rwr.normalize.affinity logical whether normalize the activity (affinity) result score using quantile normalisation,
#' default is FALSE.
#' @param rwr.threads the threads to run Random Walk With Restart (RWR), default is 2L.
#' @param hyper.test.weighted logical whether consider weighting the enrichment score of cell using hypergeometric test,
#' default is TRUE.
#' @param sv.used.reduction character used as spatial coordinates to detect SVG, default is \code{UMAP},
#' if \code{data} has \code{spatialCoords}, which will be used as spatial coordinates.
#' @param sv.grid.n numeric number of grid points in the two directions to estimate 2D weighted kernel density, default is 100.
#' @param sv.permutation numeric the number of permutation for each single feature to detect the signicantly spatially or single 
#' cell variable features, default is 100.
#' @param sv.p.adjust.method character the method to adjust the pvalue of the result, default is \code{BY}.
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
#' @return a \linkS4class{SVPExperiment} or a \linkS4class{SingleCellExperiment}, see details.
#' @details
#' if input is a \linkS4class{SVPExperiment}, output will be also a \linkS4class{SVPExperiment}, the activity score of gene sets
#' was stored in \code{assay} slot of the specified \code{gsvaexp}, and the spatially variable gene sets result is stored in \code{svDfs} 
#' of the specified \code{gsvaexp}, which is a \linkS4class{SingleCellExperiment}. If input is a \linkS4class{SingleCellExperiment} 
#' (which is extracted from \linkS4class{SVPExperiment} using \code{gsvaExp()} funtion), output will be also a
#' \linkS4class{SingleCellExperiment}, the activity score of gene sets result can be extracted using \code{assay()} function, and the
#' spatial variable gene sets result can be extracted using \code{svDf} function.
#' The result of \code{svDf} will return a matrix which has \code{sp.kld}, \code{boot.sp.kld.mean}, \code{boot.sp.kld.sd}, \code{pvalue},
#' \code{padj} and \code{rank}.
#' \itemize{
#'   \item \code{sp.kld} which is logarithms of Kullback–Leibler divergence, larger value meaning the greater the difference from the 
#'      background distribution without spatial variability.
#'   \item \code{boot.sp.kld.mean} which is mean of logarithms of Kullback–Leibler divergence based on the permutation of each features.
#'   \item \code{boot.sp.kld.sd} which is standard deviation of logarithms of Kullback–Leibler divergence based on the permutation of 
#'      each features.
#'   \item \code{pvalue} the pvalue is calculated using the real \code{sp.kld} and the permutation \code{boot.sp.kld.mean} and 
#'      \code{boot.sp.kld.sd} based on the normal distribution.
#'   \item \code{padj} the adjusted pvalue based on the speficied \code{sv.p.adjust.method}, default is \code{BY}.
#'   \item \code{rank} the order of significant spatial variable features based on \code{padj} and \code{sp.kld}.
#' }
#' @seealso [`sc.rwr`] to calculate the activity score of gene sets and [`kldSVG`] to identify the spatiall variable or specified 
#' cell gene sets or a features.
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
    knn.consider.spcoord = FALSE,
    sp.alpha.add.weight = .2,
    sp.beta.add.mp.weight = .1,    
    knn.used.reduction.dims = 30,
    knn.combined.cell.feature = FALSE,
    knn.graph.weighted = TRUE,
    knn.k.use = 400,
    rwr.restart = .75,
    rwr.normalize.adj.method = c("laplacian", "row", "column", "none"),
    rwr.normalize.affinity = FALSE,
    rwr.threads = 2L,
    hyper.test.weighted = TRUE,
    sv.used.reduction = c('UMAP', 'TSNE'),
    sv.grid.n = 100,
    sv.permutation = 100,
    sv.p.adjust.method = "BY",
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
    knn.consider.spcoord = FALSE,
    sp.alpha.add.weight = .2,
    sp.beta.add.mp.weight = .1,    
    knn.used.reduction.dims = 30,
    knn.combined.cell.feature = FALSE,    
    knn.graph.weighted = TRUE,
    knn.k.use = 400,
    rwr.restart = .75,
    rwr.normalize.adj.method = c("laplacian", "row", "column", "none"),
    rwr.normalize.affinity = FALSE,
    rwr.threads = 2L,
    hyper.test.weighted = TRUE,
    sv.used.reduction = c('UMAP', 'TSNE'),
    sv.grid.n = 100,
    sv.permutation = 100,
    sv.p.adjust.method = "BY",
    sv.BPPARAM = SerialParam(),
    run.sv = TRUE,
    cells = NULL,
    features = NULL,
    verbose = TRUE,
    random.seed = 1024,
    ...
  ){
    sce <- sc.rwr(data, 
              gset.idx.list, 
              gsvaExp.name, 
              min.sz, 
              max.sz, 
              gene.occurrence.rate, 
              assay.type, 
              knn.consider.spcoord, 
              sp.alpha.add.weight,
              sp.beta.add.mp.weight,
              knn.used.reduction.dims,
              knn.combined.cell.feature,
              knn.graph.weighted,
              knn.k.use,
              rwr.restart,
              rwr.normalize.adj.method,
              rwr.normalize.affinity,
              rwr.threads,
              hyper.test.weighted,
              cells,
              features,
              verbose
        )
    if (run.sv){
        sce <- kldSVG(data = sce, 
                  sv.used.reduction = sv.used.reduction,
                  sv.grid.n = sv.grid.n,
                  sv.permutation = sv.permutation,
                  sv.p.adjust.method = sv.p.adjust.method,
                  verbose = verbose,
                  random.seed = random.seed,
                  gsvaexp = gsvaExp.name)
    }
    return(sce)
})
