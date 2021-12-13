kms.mod <- function (x, y, H, max.iter = 400, kernel = "Fixed", K, tol.iter, tol.clust, min.clust.size, 
    merge = TRUE, keep.path = FALSE, verbose = FALSE) 
{
    n <- nrow(x)
    d <- ncol(x)
    if (missing(tol.iter)) 
        tol.iter <- 0.001 * min(apply(x, 2, IQR, na.rm=T),na.rm=T)
    if (missing(tol.clust)) 
        tol.clust <- 0.01 * max(apply(x, 2, IQR,na.rm=T),na.rm=T)
    if (missing(y)) 
        y <- x
    if (missing(min.clust.size)) 
        min.clust.size <- round(0.01 * nrow(y), 0)
    if (missing(H)) 
        H <- Hpi(x, deriv.order = 1, binned = default.bflag(d = d, n = n), nstage = 2 - (d > 2))
    Hinv <- chol2inv(chol(H))
    if (is.vector(y)) 
        y <- matrix(y, nrow = 1)
    n.seq <- block.indices(n, nrow(y), d = d, r = 0, diff = FALSE, block.limit = 3e+06)
    if (verbose) 
        pb <- txtProgressBar()
    ms <- list()
    i <- 1
    if (verbose) 
        setTxtProgressBar(pb, i/(length(n.seq) - 1))
    #ms <- kms.base.mod(x = x, y = y[n.seq[i]:(n.seq[i + 1] - 1), ], H = H, kernel = kernel, K = K, tol.iter = tol.iter, tol.clust = tol.clust, 
    #    Hinv = Hinv, verbose = verbose, max.iter = max.iter)
    ms <- kms.base.mod(x = x, y = y, H = H, kernel = kernel, K = K, tol.iter = tol.iter, tol.clust = tol.clust, 
        Hinv = Hinv, verbose = verbose, max.iter = max.iter)
    
    if (length(n.seq) < 0) {
        for (i in 2:(length(n.seq) - 1)) {
            if (verbose) 
                setTxtProgressBar(pb, i/(length(n.seq) - 1))
            ms.temp <- kms.base.mod(x = x, y = y[n.seq[i]:(n.seq[i + 
                1] - 1), ], H = H, kernel = kernel, K = K, tol.iter = tol.iter, tol.clust = tol.clust, 
                Hinv = Hinv, verbose = verbose, max.iter = max.iter)
            ms$y <- rbind(ms$y, ms.temp$y)
            ms$end.points <- rbind(ms$end.points, ms.temp$end.points)
            ms$label <- c(ms$label, ms.temp$label + max(ms$label))
            ms$mode <- rbind(ms$mode, ms.temp$mode)
            ms$nclust <- ms$nclust + ms.temp$nclust
            ms$nclust.table <- table(ms$label)
            ms$path <- c(ms$path, ms.temp$path)
            ms <- ms.merge.dist(ms = ms, tol = tol.clust, verbose = FALSE)
        }
    }
    if (verbose) 
        close(pb)
    path.temp <- ms$path
    ms$path <- NULL
    ms$tol.iter <- tol.iter
    ms$tol.clust <- tol.clust
    ms$min.clust.size <- min.clust.size
    ms$names <- parse.name(x)
    if (keep.path) 
        ms$path <- path.temp
    if (merge) 
        ms <- ms.merge.num(ms, min.clust.size = min.clust.size, 
            verbose = verbose)
    return(ms)
}


kms.base.mod <- function (x, H, Hinv, y, kernel, K, max.iter, tol.iter, tol.clust, verbose = FALSE) 
{
    if (!is.matrix(x)) 
        x <- as.matrix(x)
    if (!is.matrix(y)) 
        y <- as.matrix(y)
    if (missing(Hinv)) 
        Hinv <- chol2inv(chol(H))

    # Compute distances
    dist_kNN <- (FNN::get.knnx(x, x, k = K)$nn.dist[, K])
    dist_kNN[dist_kNN == 0] <- min(dist_kNN[dist_kNN > 0])

    nx <- nrow(x)
    ny <- nrow(y)
    d <- ncol(y)
    y.path <- split(y, row(y), drop = FALSE)
    names(y.path) <- NULL

    xHinv <- x %*% Hinv
    xHinvx <- rowSums(xHinv * x)

    y.update <- y
    i <- 1
    eps <- max(sqrt(rowSums(y.update^2)))
    disp.ind <- head(sample(1:nrow(y)), n = min(100, nrow(y)))
    while (eps > tol.iter & i < max.iter) {
        y.curr <- y.update
        yHinvy <- t(rowSums(y.curr %*% Hinv * y.curr))
        if (kernel == "Fixed") {
          Mah <- apply(yHinvy, 2, "+", xHinvx) - 2 * xHinv %*% t(y.curr)
        } else if (kernel == "Adaptive") {
          Mah <- Mah <- (apply(yHinvy, 2, "+", xHinvx) - 2 * xHinv %*% t(y.curr))/dist_kNN
        }
        w <- exp(-Mah/2)
        denom <- colSums(w)
        num <- t(w) %*% x
        denom[denom <= 0.001 * tol.iter] <- 0.001 * tol.iter
        mean.shift.H <- num/denom
        y.update <- mean.shift.H
        y.update.list <- split(y.update, row(y.update), drop = FALSE)
        y.path <- mapply(rbind, y.path, y.update.list, SIMPLIFY = FALSE)
        eps <- max(sqrt(rowSums((y.curr - y.update)^2)))

        if (verbose > 1) {
            y.range <- apply(y, 2, range)
            if (d == 2) 
                plot(y.update[disp.ind, ], col = 1, xlim = y.range[, 
                  1], ylim = y.range[, 2], xlab = "x", ylab = "y")
            else pairs(y.update[disp.ind, ], col = 1)
        }
        i <- i + 1
    }
    ms.endpt <- t(sapply(y.path, tail, n = 1, SIMPLIFY = FALSE))
    mode.tree <- hclust(dist(ms.endpt))
    clust.label <- cutree(mode.tree, h = tol.clust)
    nclust <- length(unique(clust.label))
    mode.val <- by(ms.endpt, INDICES = clust.label, FUN = colMeans)
    mode.val <- t(sapply(mode.val, FUN = identity))
    colnames(mode.val) <- colnames(x)
    rownames(mode.val) <- NULL
    nclust.table <- table(clust.label, dnn = "")
    ms <- list(x = x, y = y, end.points = ms.endpt, H = H, label = clust.label, 
        nclust = nclust, nclust.table = nclust.table, mode = mode.val, 
        path = y.path)
    class(ms) <- "kms"
    return(ms)
}

parse.name <- function (x) 
{
    if (is.vector(x)) {
        d <- 1
        x.names <- deparse(substitute(x))
    }
    else {
        d <- ncol(x)
        x.names <- colnames(x)
        if (is.null(x.names)) {
            x.names <- strsplit(deparse(substitute(x)), "\\[")[[1]][1]
            x.names <- paste(x.names, "[, ", 1:d, "]", sep = "")
        }
    }
    return(x.names)
}

ms.merge.dist <- function (ms, tol, verbose) 
{
    if (missing(tol)) 
        tol <- 0.1 * min(apply(ms$x, 2, IQR))
    mode.tree <- hclust(dist(ms$mode))
    merge.label <- cutree(mode.tree, h = tol)
    merge.label <- split(1:ms$nclust, merge.label, drop = FALSE)
    ms.temp <- ms.merge.label(ms = ms, label = merge.label, verbose = verbose)
    return(ms.temp)
}

ms.merge.label <- function (ms, label, verbose = FALSE) 
{
    ms.merge <- ms
    for (i in 1:length(label)) {
        labeli <- label[[i]]
        if (length(labeli) > 1) {
            merge.label <- min(labeli)
            ms.merge$label[ms$label %in% labeli] <- merge.label
            mode.label <- ms$mode[labeli[round(length(labeli)/2, 0)], ]
            for (j in labeli) ms.merge$mode[j, ] <- mode.label
        }
    }
    ms.merge$mode <- unique(ms.merge$mode)
    ms.merge$nclust <- nrow(ms.merge$mode)
    ms.merge$label <- as.factor(ms.merge$label)
    levels(ms.merge$label) <- 1:ms.merge$nclust
    ms.merge$label <- as.numeric(ms.merge$label)
    ms.merge$nclust.table <- table(ms.merge$label)
    if (verbose) {
        cat("Current clusters:", ms.merge$nclust.table, "\n")
    }
    return(ms.merge)
}

block.indices <- function (nx, ny, d, r = 0, diff = FALSE, block.limit = 1e+06, 
    npergroup) 
{
    if (missing(npergroup)) {
        if (diff) 
            npergroup <- max(c(block.limit%/%(nx * d^r), 1))
        else npergroup <- max(c(block.limit%/%nx, 1))
    }
    nseq <- seq(1, ny, by = npergroup)
    if (tail(nseq, n = 1) <= ny) 
        nseq <- c(nseq, ny + 1)
    if (length(nseq) == 1) 
        nseq <- c(1, ny + 1)
    return(nseq)
}

ms.merge.num <- function (ms, min.clust.size, verbose = FALSE) 
{
    if (missing(min.clust.size)) 
        min.clust.size <- round(0.01 * nrow(ms$y), 0)
    min.clust.size <- round(min.clust.size, 0)
    if (any(ms$nclust.table <= min.clust.size)) {
        if (verbose) 
            cat("Min cluster size merging begins. Min size = ", 
                min.clust.size, "\n")
        ms.temp <- ms
        while (any(ms.temp$nclust.table <= min.clust.size) & 
            ms.temp$nclust > 1) {
            nclust.table <- table(ms.temp$label)
            small.clust.ind <- which.min(nclust.table)
            if (nclust.table[small.clust.ind] <= min.clust.size) {
                nearest.clust.ind <- FNN::get.knnx(ms.temp$mode, 
                  ms.temp$mode, k = 2)$nn.index[small.clust.ind, 2]
                merge.label <- ms$label
                merge.label[merge.label == small.clust.ind] <- nearest.clust.ind
                merge.label <- split(ms$label, merge.label, drop = FALSE)
                merge.label <- lapply(merge.label, unique)
                ms.temp <- ms.merge.label(ms = ms.temp, label = merge.label, verbose = verbose)
            }
        }
        ms <- ms.temp
        ms$min.clust.size <- min.clust.size
        if (verbose) 
            cat("Min cluster size merging ends.\n\n")
    }
    if (verbose) {
        cat("Final clusters:\n")
        summary(ms)
    }
    return(ms)
}

NNMS <- function (data, K, tol, max.iter = 500, tol.iter) {
  if (missing(tol.iter)) {
    tol.iter <- 0.001 * min(apply(data, 2, IQR,na.rm=T),na.rm=T)
  }
  
  new_values <- data
  eps <- max(sqrt(rowSums(new_values^2)),na.rm=T)
  i <- 1
  
  while (eps > tol.iter & i < max.iter) {
    old_values <- new_values
    neighbors <- knn_meanShift(new_values, data, k = K)
    new_values <- t(apply(neighbors$neighbors, 1, function(x) colMeans(data[x, ])))
    eps <- max(sqrt(rowSums((new_values - old_values)^2)))
    #print(i)
    i <- i + 1
  }
  
  new_values_dist <- dist(new_values)
  dend <- hclust(new_values_dist)
  clust <- cutree(dend, h = tol)
  return(list(group = clust, final_step = new_values))
}

try2 <- function(code, silent = FALSE) {
  tryCatch(code, error = function(c) {
    if (!silent) {res$label = NA}
    else {code}})}