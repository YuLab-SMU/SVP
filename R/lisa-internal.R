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
      nm <- c("Gi", "E.Gi", "Var.Gi", "Z.Gi")
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
        x$`cluster.no.test` <- cut(x$`Z.Gi`,
                         c(-Inf, mean(x$`Z.Gi`, na.rm=TRUE), Inf),
                         lbs)
    }else{
        c1 <- cut(x$z, c(-Inf, mean(x$z), Inf), lbs)
        c2 <- cut(x$lz, c(-Inf, mean(x$lz), Inf), lbs)
        x$`cluster.no.test` <- interaction(c1, c2, sep="-")
        x$z <- NULL
        x$lz <- NULL
    }
    x$`cluster.test`<- ifelse(x[[5]] <= 0.05, as.character(x[["cluster.no.test"]]), "NoSign")
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
    nm <- names(res) |> unique()
    res <- lapply(nm, function(x)dplyr::bind_rows(res[names(res)==x]))
    names(res) <- nm
}
