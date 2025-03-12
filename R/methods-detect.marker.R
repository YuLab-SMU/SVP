#' @title Detecting the specific cell features with nearest distance of cells in MCA space
#' @rdname runDetectMarker-method
#' @param data SingleCellExperiment object
#' @param group.by the column name of cell annotation. Or a vector of length equal to 
#' \code{ncol(data)}, specifying the group to which each cell is assigned. It is required.
#' @param aggregate.group logical whether calculate the center cluster of each group of cell according
#' to the \code{group.by}, then find the nearest features of the center cluster, default TRUE. If
#' FALSE, meaning the nearest features to each cell are detected firstly.
#' @param reduction character which reduction space, default is \code{'MCA'}.
#' @param dims integer the number of components to defined the nearest distance.
#' @param ntop integer the top number of nearest or furthest (\code{type = 'negative'}) features, 
#' default is 200.
#' @param present.prop.in.group numeric the appearance proportion of groups which have the marker
#' default is .1, smaller value represent the marker will have higher specificity, but the number of 
#' marker for each group might also decrease, the minimum value is \code{1/length(unique(data[[group.by]]))}.
#' @param present.prop.in.sample numeric the appearance proportion of samples which have the marker in 
#' the corresponding group by specific \code{group.by}, default is 0.2.
#' @param BPPARAM A BiocParallelParam object specifying whether perform the analysis in parallel using
#' \code{BiocParallel} default is \code{SerialParam()}, meaning no parallel.
#' You can use \code{BiocParallel::MulticoreParam(workers=4, progressbar=TRUE)} to parallel it,
#' the \code{workers} of \code{MulticoreParam} is the number of cores used, see also
#' \code{\link[BiocParallel]{MulticoreParam}}. default is \code{SerialParam()}.
#' @param ... additional parameters.
#' @return a list, which contains features and named with clusters of \code{group.by}.
#' @export
#' @examples
#' # The example data (small.sce) is generated through simulation and has no actual meaning.
#' set.seed(123)
#' example(runMCA, echo = FALSE)
#' small.sce |> runDetectMarker(group.by = 'Cell_Cycle', ntop = 20, 
#'               present.prop.in.sample = .2)
#' # group.by, a vector of length equal to ncol(small.sce)
#' small.sce |> runDetectMarker(
#'                group.by = small.sce$Cell_Cycle, 
#'                ntop = 20,
#'                present.prop.in.sample = .2
#'              )
setGeneric('runDetectMarker',
  function(
    data,
    group.by,
    aggregate.group = TRUE,
    reduction = 'MCA',
    dims = 30,
    ntop = 200,
    present.prop.in.group = 0.1,
    present.prop.in.sample = .2,
    BPPARAM = SerialParam(),
    ...
  )
  standardGeneric('runDetectMarker')
)

#' @rdname runDetectMarker-method
#' @aliases runDetectMarker,SingleCellExperiment
#' @export runDetectMarker
#' @importFrom stats setNames
setMethod(
  'runDetectMarker', 
  'SingleCellExperiment', 
  function(
    data, 
    group.by,
    aggregate.group = TRUE,
    reduction = 'MCA', 
    dims = 30, 
    ntop = 200,
    present.prop.in.group = .1,
    present.prop.in.sample = .2,
    BPPARAM = SerialParam(),
    ...
  ){
  
    rd <- reducedDim(data, reduction)
    if (grepl("MCA", reduction, ignore.case=TRUE)){
        rd.f.nm <- "genesCoordinates"
    }
    if (grepl("PCA", reduction, ignore.case=TRUE)){
        rd.f.nm <- "rotation"
    }

    f.rd <- attr(rd, rd.f.nm)
    if(length(dims)==1){
        dims <- seq(min(ncol(rd), dims))
    }else{
        dims <- seq(min(ncol(rd), max(dims)))
    }
    
    cell.rd <- rd[, dims, drop=FALSE]

    f.rd <- f.rd[,dims, drop=FALSE]

    if (.check_group.by(group.by, ncol(data))){
        group.vec <- as.character(group.by)
    }else{
        group.vec <- as.character(colData(data)[[group.by]])
    }

    if (aggregate.group){
        cell.rd <- .calGroupCenter(cell.rd, group.vec, fun = "mean", BPPARAM)
    }

    
    dt <- .build_pair_knn(cell.rd, f.rd, ntop, FALSE)

    if (present.prop.in.sample > 1){
        present.prop.in.sample <- .5
    }
    
    gsetlist <- .check.genes.present2(
                   data,
                   dt, 
                   group.vec, 
                   present.prop.in.group,
                   present.prop.in.sample,
                   BPPARAM
                )
    return(gsetlist)
})

.calGroupCenter <- function(da, group.vec, fun="mean", BPPARAM = SerialParam()){
    fun <- switch(fun, mean=colMeans, sum=colSums)
    nm <- unique(group.vec)
    da <- bplapply(nm, function(x){
        tmpda <- da[group.vec == x, ,drop=FALSE]
        tmpda <- fun(tmpda)
        return(tmpda)
    }, BPPARAM = BPPARAM) 
    da <- do.call('rbind', da)
    rownames(da) <- nm
    return(da)
}

#.obtain.nn.genes <- function(d, type = c('positive'), ntop = 200){
#    if (type == 'positive'){
#        dt <- .build.adj.m(d, top.n = ntop)
#        return(dt)
#    }
#    d2 <- .normalize_dist_as_matrix(d)
#    if (type == 'negative'){
#        dt <- .build.adj.m(d2, top.n = ntop)
#    }else{
#        dt1 <- .build.adj.m(d, top.n = ntop)
#        dt2 <- .build.adj.m(d2, top.n = ntop)
#        dt <- dt1 + dt2
#    }
#    return(dt)
#}

.normalize_dist_as_matrix <- function(x){
    y <- .normalize_dist(x)
    attr(y, "dim") <- attr(x, "dim")
    attr(y, "dimnames") <- attr(x, "dimnames")
    return(y)
}

.check.genes.present <- function(da, group.vec, pseudomarker, present.prop = .5, BPPARAM){
    nm <- names(pseudomarker)
    da <- assay(da)
    
    res <- bplapply(nm, function(i){
       tmpda <- da[pseudomarker[[i]], group.vec == i, drop=FALSE]
       flag <- DelayedMatrixStats::rowMeans2(tmpda > 0) >= present.prop
       rownames(tmpda)[flag]
    }, BPPARAM = BPPARAM) |> setNames(nm)
    
    return(res)

}

.check.genes.present2 <- function(
    data,
    da, 
    group.vec, 
    present.prop.in.group, 
    present.prop.in.sample, 
    BPPARAM
    ){
    flag <- ncol(da) == length(group.vec)
    if (flag){
        da <- Matrix::t(da) |> as.matrix() |> as.data.frame(check.names=FALSE)
        da <- .calGroupCenter(da, group.vec, fun = "mean", BPPARAM)
        da <- t(da)
    }
    
    present.prop.in.group <- max(present.prop.in.group, 1/ncol(da))
    index <- DelayedMatrixStats::rowMeans2(da > 0) <= present.prop.in.group
    da <- da[index,,drop=FALSE]
    res <- apply(da, 2, function(x)names(which(x>0)), simplify=FALSE)
    names(res) <- colnames(da)
    res <- .check.genes.present(data, group.vec, res, present.prop.in.sample, BPPARAM)
    return(res)
}
