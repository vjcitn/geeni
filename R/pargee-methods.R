### =========================================================================
### pargee methods
### -------------------------------------------------------------------------
###

getDep <- function (x) deparse(formula(x)[[2]])

getEta <- function (gd, i, beta) getX(gd, i) %*% beta

getMu <- function (gd, i, beta, family) 
    as.numeric(family()$linkinv(getEta(gd, i, beta)))

getX <- function (gd, i) 
{
    dat <- getGrp(gd, i)
    dep <- getDep(gd)
    cbind(1,dat[, colnames(dat) != dep, drop = FALSE])
}

getY <- function (gd, i) getGrp(gd, i)[, getDep(gd)]

Di = function (gd, i, beta, family) 
    as.numeric(family()$mu.eta(getEta(gd, i, beta))) * getX(gd, i)

Vinv.i <- function(gd, i, beta, family) 
    diag(1/family()$variance(getMu(gd, i, beta, family)))

ri <- function(gd, i, beta, family) getY(gd,i) - getMu(gd, i, beta, family)

Gcomps <- function(i, gd, beta, family) {
    DD <- Di(gd, i, beta, family)
    val <- t(DD) %*% Vinv.i(gd, i, beta, family) 
    val1 <- val %*% DD
    val2 <- val %*% ri(gd, i, beta, family)
    list(DtVDi=val1, DtVri=val2)
}

combi <- function(x) {
    rr <- list(0, 0)
    for (i in seq_along(x))
        rr <- Map("+", rr, x[[i]])
    rr
}

pargee  <- function(gd, family, binit, maxit=20, tol=1e-6, ...) {
    beta <- binit
    del <- Inf
    curit <- 1
    nclus <- length(discrim(gd))
    while (max(abs(del)) > tol ) {
        res <- bplapply(1:nclus, Gcomps, gd=gd, beta=beta, family=family)
        delcomp <- combi(res) 
        del <- solve(delcomp[[1]]) %*% delcomp[[2]]
        beta <- beta + del
        curit <- curit + 1
        if (curit > maxit) 
            stop(paste("maxit [", maxit, "] iterations exceeded"))
    }
    beta
}
