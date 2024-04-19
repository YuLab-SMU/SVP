#' @title Detecting the specific cell features with nearest distance of cells in MCA space
#' @rdname runDetectMarker-method
#' @param data SingleCellExperiment object
#' @param group.by the column name of cell annotation.
#' @param aggregate.group logical whether calculate the center cluster of each group of cell according
#' to the \code{group.by}, then find the nearest features of the center cluster, default TRUE. If
#' FALSE, meaning the nearest features to each cell are detected firstly.
#' @param reduction character which reduction space, default is \code{'MCA'}.
#' @param dims integer the number of components to defined the nearest distance.
#' @param ntop integer the top number of nearest or furthest (\code{type = 'negative'}) features, 
#' default is 200 .
#' @param present.prop numeric the appearance proportion of samples in the corresponding group in \code{group.by},
#' default is 0.2 .
#' @param type character which features are be extracted, the nearest features (\code{type='positive'}) 
#' or furthest features (\code{type = 'negative'}) or both (\code{type='all'}), default is 
#' \code{type='positive'}.
#' @param consider.unique.in.group logical whether detect the unique features belonging to one cell cluster.
#' @param ... additional parameters.
#' @return a list, which contains features and named with clusters of \code{group.by}.
#' @export
#' @examples
#' example(runMCA, echo = FALSE)
#' small.sce |> runDetectMarker(group.by='Cell_Cycle', ntop=20, 
#'               present.prop=.9, consider.unique.in.group=FALSE)
setGeneric('runDetectMarker',
  function(
    data,
    group.by,
    aggregate.group = TRUE,
    reduction = 'MCA',
    dims = 30,
    ntop = 200,
    present.prop = .2,
    type = c('positive', 'all', 'negative'),
    consider.unique.in.group = TRUE,
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
    present.prop = .2,
    type = c('positive', 'all', 'negative'),
    consider.unique.in.group = TRUE,
    ...
  ){
  
    type <- match.arg(type)
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
    group.vec <- as.character(colData(data)[[group.by]])

    if (aggregate.group){
        cell.rd <- as.data.frame(cell.rd, check.names=FALSE)
        cell.rd$Group <- group.vec
        cell.rd <- .calGroupCenter(cell.rd, group.by = "Group", fun = mean)
    }

    cell2features.dist <- t(pairDist(f.rd, cell.rd))
    
    dt <- .obtain.nn.genes(cell2features.dist, type, ntop)
    if (!aggregate.group){
        gnm <- unique(group.vec)
        dt <- lapply(gnm, function(x) unique(c(dt[, x==group.vec]))) |> setNames(gnm)
    } 
    if (consider.unique.in.group){
        if (aggregate.group){
            gsetlist <- lapply(seq(ncol(dt)),function(i)setdiff(dt[,i], c(dt[,-i]))) |>
                setNames(colnames(dt))
        }else{
            gsetlist <- lapply(seq(length(dt)), function(i)setdiff(dt[[i]], unlist(dt[-i]))) |>
                setNames(names(dt))
        }
    }else{
        if (aggregate.group){
            gsetlist <- lapply(seq(ncol(dt)), function(i)dt[,i]) |> setNames(colnames(dt))
        }else{
            gsetlist <- dt
        }
    }
    
    if (present.prop > 1){
        present.prop <- .5
    }
    gsetlist <- .check.genes.present(data, group.by, gsetlist, present.prop)
    return(gsetlist)
})

.calGroupCenter <- function(da, group.by, fun=mean){
    da <- da |> dplyr::group_by(!!as.symbol(group.by)) |> 
        dplyr::summarize_all(fun, na.rm=TRUE) |> 
        as.data.frame(check.names=FALSE)
    rownames(da) <- da[[group.by]]
    da[[group.by]] <- NULL
    return(as.matrix(da))
}

.obtain.nn.genes <- function(d, type = c('positive'), ntop = 200){
    if (type == 'positive'){
        dt <- apply(d, 2, function(x)names(sort(x, method = 'quick'))[seq(ntop)])
    }else if (type == 'negative'){
        dt <- apply(d, 2, function(x)names(sort(x, method = 'quick', decreasing=TRUE))[seq(ntop)])
    }else if (type == 'all'){
        dt1 <- apply(d, 2, function(x)names(sort(x, method = 'quick'))[seq(ntop)])
        dt <- apply(d, 2, function(x)names(sort(x, method = 'quick', decreasing=TRUE))[seq(ntop)])
        dt <- rbind(dt1, dt)
    }
    return(dt)
}



.check.genes.present <- function(da, group.by, pseudomarker, present.prop = .5){
    nm <- names(pseudomarker)
    
    res <- lapply(nm, function(i){
       tmpda <- da[pseudomarker[[i]], da[[group.by]]==i,drop=FALSE]
       flag <- Matrix::rowMeans(assay(tmpda) > 0) >= present.prop
       rownames(tmpda)[flag]
    }) |> setNames(nm)
    
    return(res)

}
