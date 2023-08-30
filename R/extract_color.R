.extract_color <- function(img, beta = 2, coords){
    img.rgb <- .convert.rgb(img)
    step <- round(beta/2)
    coords <- round(coords)
    img.rgb.m <- .extract.region(img.rgb, coords, step, img)
    res <- .cal_weights_rgb_spots(img.rgb.m)
    return(res)
}

.cal_weights_rgb_spots <- function(x){
    y <- matrix(apply(x, 1, var), nrow = 1)
    y <- t(y %*% x) / sum(y)
    colnames(y) <- 'weights.rbg.spots'
    return(y)
}

.convert.rgb <- function(x){
    x <- apply(x, 2, col2rgb, simplify = FALSE)
    r.x <- do.call('rbind', lapply(x, function(i)i[1, ])) |> t()
    g.x <- do.call('rbind', lapply(x, function(i)i[2, ])) |> t()
    b.x <- do.call('rbind', lapply(x, function(i)i[3, ])) |> t()
    list(r.x, g.x, b.x)
}

.extract.region <- function(x, coords, step, img){
    min.max <- .generate_min_max(coords, step, img)
    mapply(.cal_mean_color, 
           apply(min.max[[1]],1,
                 function(x)x,simplify=FALSE), 
           apply(min.max[[2]],1,                                  
                function(x)x,simplify = FALSE), 
           MoreArgs = list(x = x)
    )
}

.cal_mean_color <- function(x, range.x, range.y, fun.internal = mean){
    do.call('rbind', lapply(x, function(i)fun.internal(i[seq(range.x[2], range.y[2]), seq(range.x[1], range.y[1])])))
}


.generate_min_max <- function(x, step, img){
    if (inherits(img, 'raster')){
        width <- dim(img)[2]
        height <- dim(img)[1]
    }else if (inherits(img, 'numeric')){
        width <- img[2]
        height <- img[1]
    }
    min.x <- x - step
    max.x <- x + step
    min.x[min.x < 0] <- 0
    max.x[max.x[,1] > width, 1] <- width
    max.x[max.x[,2] > height, 2] <- height
    return(list(min.x, max.x))
}


