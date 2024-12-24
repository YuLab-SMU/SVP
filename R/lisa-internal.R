#' @importFrom stats pnorm
.internal.runLISA <- function(
    x,
    weight,
    method = c("localG", "localmoran"),
    flag.method = c("mean", "median"),
    alternative = c("two.sided", "greater", "less"),
    BPPARAM = SerialParam()
  ){
  if (!inherits(x, 'dgCMatrix')){
    x <- as.matrix(x) |> Matrix::Matrix(sparse = TRUE) 
  }  
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  flag.method <- match.arg(flag.method)
  prefix <- switch(alternative, two.sided = "!=", greater = ">", less = "<")

  if (method == 'localG'){
      nm <- c("Gi", "E.Gi", "Var.Gi", "Z.Gi", "x")
      pnm <- gettextf("Pr (z %s E(%s))", prefix, "Gi")
      nm <- c(nm, pnm)
      res <- CalLocalGParallel(x, weight) 
  }else{
      nm <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", "z", "lz")
      pnm <- gettextf("Pr (z %s E(%s))", prefix, "Ii")
      nm <- c(nm, pnm)
      res <- CalLocalMoranParallel(x, weight) 
  }

  res <- .add_localisa_pvalue(res, pnm, nm, colnames(x), method, flag.method, alternative, BPPARAM)
  names(res) <- rownames(x)
  return(res)
}

.add_localisa_pvalue <- function(res, pnm, nm, cellnm, method, flag.method, alternative, BPPARAM){
  res <- BiocParallel::bplapply(res, function(x){
            x <- data.frame(x)
            if (alternative == 'two.sided'){
                x[[pnm]] <- 2 * pnorm(abs(x[[4]]), lower.tail = FALSE)
            }else if (alternative == 'greater'){
                x[[pnm]] <- pnorm(x[[4]], lower.tail = FALSE)
            }else{
                x[[pnm]] <- pnorm(x[[4]])
            }
            colnames(x) <- nm
            rownames(x) <- cellnm
            x <- .add_cluster(x, method, flag.method)
            return(x)
         }, BPPARAM = BPPARAM)
  return(res)
}

.add_cluster <- function(x, method, flag.method = 'mean'){
    lbs <- c("Low", "High")
    if (method == 'localG'){
        if (!all(is.na(x$`Z.Gi`))){
            x$`cluster.no.test` <- cut(x$`Z.Gi`,
                             c(-Inf, do.call(flag.method,list(x$Gi, na.rm=TRUE)), Inf),
                             lbs)
            x$x <- NULL
        }
    }else{
        if (!all(is.na(x$z)) && !all(is.na(x$lz))){
            c1 <- cut(x$z, c(-Inf, do.call(flag.method, list(x$z, na.rm=TRUE)), Inf), lbs)
            c2 <- cut(x$lz, c(-Inf, do.call(flag.method, list(x$lz, na.rm=TRUE)), Inf), lbs)
            x$`cluster.no.test` <- interaction(c1, c2, sep="-")
        }
        x$z <- NULL
        x$lz <- NULL
    }
    if (!all(is.na(x[[5]]))){
        x$`cluster.test`<- ifelse(x[[5]] <= 0.05, as.character(x[["cluster.no.test"]]), "NoSign")
    }
    return(x)
}

.add_LISA_res <- function(x, lisa.res, method){
    nm <- paste0(method,".SVP")
    lisa.res <- lapply(lisa.res, function(x)DataFrame(I(x))) |> 
                  DataFrame() |> 
                  setNames(names(lisa.res))
    if (is.null(x@int_colData$localResults)){
        x@int_colData$localResults <- matrix(nrow=ncol(x), ncol=0) |> DataFrame()
    }else{
        if (!is.null(x@int_colData$localResults[[nm]])){
            lisa.res <- .update_lisa(lisa.res, x@int_colData$localResults[[nm]])
        }
    }
    x@int_colData$localResults[[nm]] <- lisa.res
    return(x)
}

.tidy_lisa_res <- function(res){
    if (length(res)==1){
        return(res[[1]])
    }

    res <- unlist(unname(res), recursive=FALSE)
    if (length(res) == 1){
        return(res)
    }
    nm <- names(res) |> unique()
    res <- lapply(nm, function(x)dplyr::bind_rows(res[names(res)==x]))
    names(res) <- nm
    return(res)
}

.internal.as_tbl_df <- function(x, 
                                diag = FALSE, 
                                rmrd = TRUE,
                                flag.clust = FALSE, 
                                dist.method = 'euclidean', 
                                hclust.method = 'average'
    ){
    if (identical(rownames(x), colnames(x)) && nrow(x)==ncol(x) && rmrd){
        if (flag.clust && nrow(x) > 2){
            x <- .adjust_order_by_clust1(x, dist.method, hclust.method)
        }
        x <- .as_from_square_matrix(x = x, diag = diag)
        return(x)
    }
    if (flag.clust && nrow(x) > 2 && ncol(x) > 2){
        x <- .adjust_order_by_clust2(x, dist.method, hclust.method)
    }
    x <- x |>
         as.data.frame(check.names=FALSE) |>
         tibble::rownames_to_column(var='.Term') |>
         tidyr::pivot_longer(cols=-".Term") |>
         dplyr::mutate(.Term=factor(.data$.Term, rownames(x)), name=factor(.data$name, colnames(x)))
    return(x)
}

.as_from_square_matrix <- function(x, clnm = NULL, diag = FALSE){
    ind <- which(upper.tri(x, diag = diag), arr.ind = TRUE)
    nn <- rownames(x)
    x <- data.frame(x = nn[ind[,1]] , y= nn[ind[,2]], val=x[ind])
    x$x <- factor(x$x, nn)
    x$y <- factor(x$y, nn)
    if (!is.null(clnm)){
        colnames(x) <- clnm
    }
    return(x)
}

.adjust_order_by_clust1 <- function(x, dist.method, hclust.method){
    res <- dist(x, dist.method) |> hclust(hclust.method)
    x <- x[res$order,res$order,drop=FALSE]
    return(x)
}

.adjust_order_by_clust2 <- function(x, dist.method, hclust.method){
    res1 <- dist(x, dist.method) |> hclust(hclust.method)
    res2 <- dist(t(x), dist.method) |> hclust(hclust.method)
    x <- x[res1$order, res2$order, drop=FALSE]
    return(x)
}

.as_list_to_tbl <- function(x){
    f <- unlist(x)
    term <- lapply(names(x), function(i)rep(i, length(x[[i]]))) |> unlist()
    data.frame(f=f, term=term)
}

.add_list <- function(x, y){
    y <- .as_list_to_tbl(y)
    ind1 <- "f" |> setNames("x")
    ind2 <- "f" |> setNames("y")
    x <- dplyr::left_join(x, y, by = ind1)
    x <- dplyr::left_join(x, y, by = ind2)
    return(x)
}

.update_lisa <- function(x, y){
    oldn <- names(y)
    newn <- names(x)
    return(cbind(y[,setdiff(oldn, newn),drop=FALSE], x))
}
