### =========================================================================
### pargee methods
### -------------------------------------------------------------------------
###

### probably stay in R
getDep <- function (x) deparse(formula(x)[[2]])

getMu <- function (gd, i, beta, family) 
    as.numeric(family()$linkinv(getEta(gd, i, beta)))

getY <- function (gd, i) getGroup(gd, i)[, getDep(gd)]

getX <- function (gd, i) 
{
    dat <- getGroup(gd, i)
    dep <- getDep(gd)
    X <- cbind(1, dat[, colnames(dat) != dep, drop = FALSE])
    if (!is(gd, "gdManager"))
        as.matrix(X)
    else
        X
}

ri <- function(gd, i, beta, family) getY(gd,i) - getMu(gd, i, beta, family)

### possibly move to C
getEta <- function (gd, i, beta) getX(gd, i) %*% beta

Di <- function (gd, i, beta, family) 
    as.numeric(family()$mu.eta(getEta(gd, i, beta))) * getX(gd, i)

Vinv.i <- function(gd, i, beta, family) 
    diag(1/family()$variance(getMu(gd, i, beta, family)))

Gcomps <- function(i, gd, beta, family, sandwich=TRUE) {
    DD <- Di(gd, i, beta, family)
    tDD <- t(DD)
    r_i <- ri(gd, i, beta, family)
    Vinv <- Vinv.i(gd, i, beta, family)
    val <- tDD %*% Vinv 
    val1 <- val %*% DD
    val2 <- val %*% r_i
    middle <- NA 
    if (sandwich) {
        m1 <- (tDD %*% Vinv) %*% r_i
        m2 <- r_i %*% (Vinv %*% DD) 
        middle <- m1 %*% m2
    }
    list(DtVDi=val1, DtVri=val2, DtVririVD=middle)
}

combi <- function(x) {
    xx <- x[[1]] 
    for (i in 2:length(x))
        xx <- Map("+", xx, x[[i]]) 
    xx
}

setGeneric("pargee", 
    function(gd, family, binit, maxit=20, tol=1e-6, sandwich=TRUE, ...)
       standardGeneric("pargee")
)

.pargee <- function(gd, family, binit, maxit, tol, sandwich, ngroup) {
    beta <- binit
    del <- Inf
    curit <- 1
    robvar <- NA
    while (max(abs(del)) > tol ) {
        res <- bplapply(1:ngroup, Gcomps, gd=gd, beta=beta, family=family, 
                        sandwich=sandwich)
        delcomp <- combi(res) 
        solve_DtVDi <- solve(delcomp[[1]])
        del <- solve_DtVDi %*% delcomp[[2]]
        beta <- beta + del
        if (sandwich)
            robvar <- solve_DtVDi %*% (delcomp[[3]] %*% solve_DtVDi) 
        curit <- curit + 1
        if (curit > maxit) 
            stop(paste("maxit [", maxit, "] iterations exceeded"))
    }
    list(coefficients=beta, robust.variance=robvar)
}
setMethod("pargee", "gdManager", 
    function(gd, family, binit, maxit=20, tol=1e-6, sandwich=TRUE, ...)
        .pargee(gd, family, binit, maxit, tol, sandwich, length(discrim(gd)))
)

setMethod("pargee", "ANY",
    function(gd, family, binit, maxit=20, tol=1e-6, sandwich=TRUE, ...)
        .pargee(gd, family, binit, maxit, tol, sandwich, 
                length(unique(getGroups(gd))))
)
