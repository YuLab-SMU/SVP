#' @importFrom stats pnorm
.internal.runLISA <- function(
    x,
    weight,
    wi,
    wi2,
    n,
    method=c("localG", "localmoran"),
    alternative = c("two.sided", "greater", "less")
  ){
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  prefix <- switch(alternative, two.sided = "!=", greater = ">", less = "<")

  if (method == 'localG'){
      nm <- c("Gi", "E.Gi", "Var.Gi", "Z.Gi", "x")
      pnm <- gettextf("Pr (z %s E(%s))", prefix, "Gi")
      nm <- c(nm, pnm)
      res <- CalLocalGCpp(x, weight, wi, wi2, n) |> data.frame()
  }else{
      nm <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", "z", "lz")
      pnm <- gettextf("Pr (z %s E(%s))", prefix, "Ii")
      nm <- c(nm, pnm)
      res <- CalLocalMoranCpp(x, weight, wi, wi2, n) |> data.frame()
  }

  if (alternative == 'two.sided'){
      res[[pnm]] <- 2 * pnorm(abs(res[[4]]), lower.tail = FALSE)
  }else if (alternative == 'greater'){
      res[[pnm]] <- pnorm(res[[4]], lower.tail = FALSE)
  }else{
      res[[pnm]] <- pnorm(res[[4]])
  }
  colnames(res) <- nm
  rownames(res) <- names(x)
  res <- .add_cluster(res, method)

  return(res)
}

.add_cluster <- function(x, method){
    lbs <- c("Low", "High")
    if (method == 'localG'){
        if (!all(is.na(x$`Z.Gi`))){
            x$`cluster.no.test` <- cut(x$`Z.Gi`,
                             c(-Inf, mean(x$x, na.rm=TRUE), Inf),
                             lbs)
            x$x <- NULL
        }
    }else{
        if (!all(is.na(x$z)) && !all(is.na(x$lz))){
            c1 <- cut(x$z, c(-Inf, mean(x$z), Inf), lbs)
            c2 <- cut(x$lz, c(-Inf, mean(x$lz), Inf), lbs)
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
    lisa.res <- lapply(lisa.res, function(x)DataFrame(I(x))) |> DataFrame() |> setNames(names(lisa.res))
    x@int_colData$localResults <- matrix(nrow=ncol(x), ncol=0) |> DataFrame()
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
