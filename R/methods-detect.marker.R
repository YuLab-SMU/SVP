setGeneric('detect.marker',
  function(
    data,
    by.group,
    reduction = 'MCA',
    dims = 30,
    ntop = 200,
    type = c('positive', 'all', 'negative'),
    consider.unique.in.group = TRUE,
    ...
  )
  standardGeneric('detect.marker')
)

#' @importFrom stats mean
setMethod(
  'detect.marker', 
  'SingleCellExperiment', 
  function(
    data, 
    by.group,
    reduction = 'MCA', 
    dims = 30, 
    ntop = 200, 
    type = c('positive', 'all', 'negative'),
    consider.unique.in.group = TRUE,
    ...
  ){
  
    type <- match.arg(type)
    rd <- reducedDim(data, reduction)
    rd.f.nm <- switch(knn.used.reduction, MCA='genesCoordinates', PCA='rotation')
    f.rd <- attr(rd, rd.f.nm)
    dims <- ifelse(length(dims)==1, seq(min(ncol(rd), dims)), seq(min(ncol(rd), max(dims))))
    
    cell.rd <- rd[, dims, drop=FALSE]
    f.rd <- f.rd[,dims, drop=FALSE]
    cell.rd$Group <- as.character(colData(data)[[by.group]])
    
    group.rd <- .calGroupCenter(cell.rd, by.group = by.group, fun = mean)

    group2features.dist <- t(pairDist(f.rd, group.rd))
    
    dt <- .obtain.nn.genes(group2features.dist, type, ntop) 
    if (consider.unique.by.group){
        gsetlist <- lapply(seq(ncol(dt)),function(i)setdiff(dt[,i], c(dt[,-i]))) |>
            setNames(colnames(dt))
    }
    return(gsetlist)
})

.calGroupCenter <- function(da, by.group, fun=mean){
    da <- da |> dplyr::group_by(!!as.symbol(by.group)) |> 
        dplyr::summarize_all(fun, na.rm=TRUE) |> 
        as.data.frame(check.names=FALSE)
    rownames(da) <- da[[by.group]]
    da[[by.group]] <- NULL
    return(da)
}

.obtain.nn.genes <- function(d, type = c('positive'), ntop = 200){
    if (type == 'positive'){
        dt <- apply(d, 2, function(x)names(SortNv(x))[seq(ntop)])
    }else if (type == 'negative'){
        dt <- apply(d, 2, function(x)names(SortNv(x, decreasing=TRUE))[seq(ntop)])
    }else if (type == 'all'){
        dt1 <- apply(d, 2, function(x)names(SortNv(x))[seq(ntop)])
        dt <- apply(d, 2, function(x)names(SortNv(x, decreasing=TRUE))[seq(ntop)])
        dt <- rbind(dt1, dt)
    }
    return(dt)
}
