gram <-
function (X, trX, block = FALSE, block.size = 1000) 
{
    if (block == FALSE) {
        if (missing(trX)) {
            tX <- t(X)
        }
        else {
            tX <- trX
        }
        G <- X %*% tX
    }
    if (block == TRUE) {
        n <- dim(X)[1]
        p <- dim(X)[2]
        if (block.size == 1) {
            stop("block.size must be greater than 1")
        }
        if (block.size >= p) {
            stop("block.size must me smaller than the number of columns of X")
        }
        if (missing(trX)) {
            tX <- t(X)
        }
        else {
            tX <- trX
        }
        col1 <- seq(1, (p - block.size), by = block.size)
        col2 <- seq(block.size, p, by = block.size)
        d <- length(col2)
        if (length(col1) != d) {
            col1[d] <- (p - block.size) + 1
        }
        col2[d] <- p
        G <- list()
        for (k in 1:d) {
            G[[k]] <- matrix(NA, n, n)
            R <- col1[k]:col2[k]
            G[[k]] <- X[, R] %*% tX[R, ]
            flush.console()
            cat("\r", paste("Gram matrix block-computation:", 
                sep = " "), paste(round(k/d * 100), "%", sep = ""))
        }
        G <- Reduce("+", G)
    }
    G
}
gram.parallel <-
function (X, cl, ...) 
{
    p <- ncol(X)
    d <- 1:p
    n.cl <- length(cl)
    batches <- split(d, ceiling(seq_along(d)/(p/n.cl)))
    Xbatches <- list()
    for (j in 1:n.cl) {
        Xbatches[[j]] <- X[, batches[[j]]]
    }
    G <- Reduce("+", parLapply(cl, Xbatches, gram, ...))
    G
}
sbr <-
function (y, X, trX, G, estimator = "PM", sparsify = FALSE, sparse.control = 1, 
    p.threshold = 5000, cov.blocks = 1000, parallel = FALSE, 
    cl, L.optim = 10^-4, U.optim = 10^4) 
{
    start.time <- proc.time()
    if (is.matrix(X) == FALSE & is.list(X) == FALSE) {
        stop("X must be either a matrix (one data-source) or a list (multiple data-sources)")
    }
    if (is.matrix(G) == FALSE & is.list(G) == FALSE) {
        stop("G must be either a matrix (one data-source) or a list (multiple data-sources)")
    }
    if (is.matrix(G) == TRUE) {
        K <- 1
    }
    if (is.list(G) == TRUE) {
        K <- length(G)
    }
    n <- length(y)
    I <- diag(n)
    ty <- t(y)
    if (missing(trX)) {
        if (K == 1) {
            tX <- t(X)
        }
        if (K != 1) {
            tX <- lapply(X, t)
        }
    }
    else {
        tX <- trX
    }
    optim.lambda <- function(args, y, G, estimator, lambda.PM.prior) {
        if (K == 1) {
            lambda <- args
            lambdaG <- lambda^(-1) * G
        }
        else {
            lambda <- c()
            lambdaG <- list()
            for (j in 1:K) {
                lambda[j] <- args[j]
                lambdaG[[j]] <- lambda[j]^(-1) * G[[j]]
            }
            lambdaG <- Reduce("+", lambdaG)
        }
        mat <- I + lambdaG
        inv.mat <- chol2inv(chol(mat))
        if (estimator == "CV") {
            opt.lambda <- log(ty %*% inv.mat %*% diag(diag(inv.mat^(-2))) %*% 
                inv.mat %*% y)
        }
        if (estimator == "ML") {
            opt.lambda <- 0.5 * determinant(mat, logarithm = TRUE)$modulus + 
                n/2 * log(0.5 * ty %*% inv.mat %*% y)
        }
        if (estimator == "PM") {
            opt.lambda <- 0.5 * determinant(mat, logarithm = TRUE)$modulus + 
                n/2 * log(0.5 * ty %*% inv.mat %*% y) + sum(lambda/lambda.PM.prior)
        }
        opt.lambda
    }
    if (estimator == "CV") {
        lambda <- optim(rep(1, K), optim.lambda, lower = rep(L.optim, 
            K), upper = rep(U.optim, K), method = "L-BFGS-B", 
            y = y, G = G, estimator = "CV")$par
    }
    if (estimator == "ML") {
        lambda <- optim(rep(1, K), optim.lambda, lower = rep(L.optim, 
            K), upper = rep(U.optim, K), method = "L-BFGS-B", 
            y = y, G = G, estimator = "ML")$par
    }
    if (estimator == "PM") {
        lambda.star <- optim(rep(1, K), optim.lambda, lower = rep(L.optim, 
            K), upper = rep(U.optim, K), method = "L-BFGS-B", 
            y = y, G = G, estimator = "CV")$par
        lambda <- optim(rep(1, K), optim.lambda, lower = rep(L.optim, 
            K), upper = rep(U.optim, K), method = "L-BFGS-B", 
            y = y, G = G, estimator = "PM", lambda.PM.prior = lambda.star)$par
    }
    inv.lambda <- lambda^(-1)
    if (K == 1) {
        GL <- inv.lambda * G
    }
    else {
        Glambda <- list()
        for (j in 1:K) {
            Glambda[[j]] <- inv.lambda[j] * G[[j]]
        }
        GL <- Reduce("+", Glambda)
    }
    IGL <- I + GL
    invIGL <- chol2inv(chol(IGL))
    yGL <- y - invIGL %*% GL %*% y
    if (K == 1) {
        beta <- inv.lambda * tX %*% yGL
    }
    else {
        beta <- list()
        for (j in 1:K) {
            beta[[j]] <- inv.lambda[j] * tX[[j]] %*% yGL
        }
        beta <- Reduce(c, beta)
    }
    RSS <- ty %*% invIGL %*% y
    sigma2 <- c(RSS/(n + 2))
    message("\n", paste("SBR Done."))
    if (sparsify == TRUE) {
        if (K == 1) {
            X <- list(X)
            tX <- list(tX)
        }
        p.source <- as.numeric(lapply(X, ncol))
        if (parallel == FALSE) {
            ind.small <- which(p.source <= p.threshold)
            ind.big <- which(p.source > p.threshold)
            var.beta <- list()
            length(var.beta) <- K
            D <- diag(1, cov.blocks)
            for (j in ind.small) {
                var.beta[[j]] <- diag(inv.lambda[j] * (diag(p.source[j]) - 
                  tX[[j]] %*% invIGL %*% X[[j]] * inv.lambda[j]))
            }
            for (j in ind.big) {
                col1 <- seq(1, (p.source[j] - cov.blocks), by = cov.blocks)
                col2 <- seq(cov.blocks, p.source[j], by = cov.blocks)
                d <- length(col2)
                if (length(col1) != d) {
                  col1[d] <- (p.source[j] - cov.blocks) + 1
                }
                col2[d] <- p.source[j]
                for (k in 1:(d - 1)) {
                  flush.console()
                  R <- col1[k]:col2[k]
                  var.beta[[j]][R] <- diag(inv.lambda[j] * (D - 
                    tX[[j]][R, ] %*% invIGL %*% X[[j]][, R] * 
                      inv.lambda[j]))
                  cat("\r", paste("Source", j, "sigma block-computation:", 
                    sep = " "), paste(round(k/(d - 1) * 100), 
                    "%", sep = ""))
                }
                k <- d
                R <- col1[k]:col2[k]
                var.beta[[j]][R] <- diag(inv.lambda[j] * (diag(1, 
                  col2[k] - col2[k - 1]) - tX[[j]][R, ] %*% invIGL %*% 
                  X[[j]][, R] * inv.lambda[j]))
                cat("\n", sep = " ")
            }
            var.beta <- Reduce(c, var.beta)
        }
        if (parallel == TRUE) {
            var.function1 <- function(X, inverse.lambda, A) {
                tX <- t(X)
                Id <- diag(ncol(X))
                variance <- diag(inverse.lambda * (Id - tX %*% 
                  A %*% X * inverse.lambda))
                variance
            }
            var.function2 <- function(X, inverse.lambda, A, n.blocks) {
                tX <- t(X)
                n.p <- ncol(X)
                variance <- c()
                if ((n.p - n.blocks) > n.blocks) {
                  Id <- diag(n.blocks)
                  col1 <- seq(1, (n.p - n.blocks), by = n.blocks)
                  col2 <- seq(n.blocks, n.p, by = n.blocks)
                  d <- length(col2)
                  if (length(col1) != d) {
                    col1[d] <- (n.p - n.blocks) + 1
                  }
                  col2[d] <- n.p
                  for (k in 1:(d - 1)) {
                    R <- col1[k]:col2[k]
                    variance[R] <- diag(inverse.lambda * (Id - 
                      tX[R, ] %*% A %*% X[, R] * inverse.lambda))
                  }
                  k <- d
                  R <- col1[k]:col2[k]
                  variance[R] <- diag(inverse.lambda * (diag(1, 
                    col2[k] - col2[k - 1]) - tX[R, ] %*% A %*% 
                    X[, R] * inverse.lambda))
                }
                else {
                  col1 <- c(1, n.blocks + 1)
                  col2 <- c(n.blocks, n.p)
                  for (k in 1:2) {
                    R <- col1[k]:col2[k]
                    variance[R] <- diag(inverse.lambda * (diag(1, 
                      length(R)) - tX[R, ] %*% A %*% X[, R] * 
                      inverse.lambda))
                  }
                }
                variance
            }
            environment(var.function1) <- .GlobalEnv
            environment(var.function2) <- .GlobalEnv
            var.beta <- list()
            n.cl <- length(cl)
            Xbatches <- list()
            n.col <- rep(NA, K)
            for (j in 1:K) {
                if (p.source[j] > n.cl) {
                  d <- 1:p.source[j]
                  batches <- split(d, ceiling(seq_along(d)/(p.source[j]/n.cl)))
                  Xbatches[[j]] <- list()
                  for (i in 1:n.cl) {
                    Xbatches[[j]][[i]] <- X[[j]][, batches[[i]]]
                  }
                  n.col[j] <- ncol(Xbatches[[j]][[1]])
                  if (n.col[j] <= p.threshold) {
                    assign(paste("Var.beta"), Reduce(c, parLapply(cl, 
                      as.array(Xbatches[[j]]), var.function1, 
                      inverse.lambda = inv.lambda[j], A = invIGL)))
                  }
                  else {
                    assign(paste("Var.beta"), Reduce(c, parLapply(cl, 
                      as.array(Xbatches[[j]]), var.function2, 
                      inverse.lambda = inv.lambda[j], A = invIGL, 
                      n.blocks = cov.blocks)))
                  }
                }
                else {
                  Var.beta <- var.function1(X[[j]], inverse.lambda = inv.lambda[j], 
                    A = invIGL)
                }
                var.beta[[j]] <- Var.beta
            }
            var.beta <- Reduce(c, var.beta)
        }
        source.weight <- lambda/sum(lambda)
        abs.beta <- abs(beta)
        kappa <- (1/abs.beta)^rep(source.weight, times = p.source)
        sparse.threshold <- RSS/n * var.beta * kappa * sparse.control
        sparse.ind <- which(abs.beta > sparse.threshold)
        sparse.beta <- rep(0, length(beta))
        sparse.beta[sparse.ind] <- beta[sparse.ind] - (beta[sparse.ind]/abs.beta[sparse.ind]) * 
            sparse.threshold[sparse.ind]
        message("\n", paste("SSBR Done."))
    }
    diff.time <- proc.time() - start.time
    setOldClass("proc_time")
    if (sparsify == FALSE) {
        setClass(Class = "Results", representation(BetaSBR = "numeric", 
            Sigma2 = "numeric", Lambda = "numeric", LambdaEstimator = "character", 
            Duration = "proc_time"))
        return(new("Results", BetaSBR = c(beta), Sigma2 = sigma2, 
            Lambda = lambda, LambdaEstimator = estimator, Duration = diff.time))
    }
    else {
        setClass(Class = "Results", representation(BetaSBR = "numeric", 
            BetaSSBR = "numeric", Sigma2 = "numeric", Lambda = "numeric", 
            LambdaEstimator = "character", Duration = "proc_time"))
        return(new("Results", BetaSBR = c(beta), BetaSSBR = sparse.beta, 
            Sigma2 = sigma2, Lambda = lambda, LambdaEstimator = estimator, 
            Duration = diff.time))
    }
}
