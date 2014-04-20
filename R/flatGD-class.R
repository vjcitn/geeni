### =========================================================================
### flatGD objects
### -------------------------------------------------------------------------
###

setOldClass("ff_vector")
setOldClass("ff_matrix")
setOldClass("formula")

setClass("gdManager", 
    representation(
        discrim = "ff_vector",
        numdat="ff_matrix", 
        formula="formula"
    )
)

### getters / setters
discrim <- function(x) x@discrim
numdat <- function(x) x@numdat
setMethod("formula", "gdManager", function(x, ...) x@formula) 

### show 
setMethod("show", "gdManager", 
    function(object) {
        cat("gdManager instance\n")
        cat("number of clusters = ", length(object@discrim), "\n")
        nr = nrow(object@numdat)
        nc = ncol(object@numdat)
        cat("size of data matrix: ", nr, "x", nc, "\n")
        cat("excerpt:\n")
        print(object@numdat[1:min(c(3,nr)),1:min(c(5,nc))])
    }
)

###  gd2flat 
gd2flat <- function(gd=Orth.new, gcol=3, prefix="gdm", overwrite=TRUE) 
{
    disc <- ff(runLength(Rle(gd[,gcol])), 
               filename=paste(prefix, "disc.ff", sep="."))
    dm <- data.matrix(gd[,-gcol])
    cn <- colnames(dm)
    rn <- rownames(gd)
    dat <- ff(data.matrix(gd[,-gcol]), dim=c(nrow(gd), ncol(gd)-1), 
              dimnames=list(rn,cn),
        filename=paste(prefix, "dat.ff", sep="."), overwrite=overwrite)
    new("gdManager", discrim=disc, numdat=dat, formula=formula(gd))
}
 
### getGrp 
setGeneric("getGrp", function(gd, ind) standardGeneric("getGrp")) 
setMethod("getGrp", c("gdManager", "numeric"), 
    function(gd, ind) {
        if (length(ind) > 1) 
            stop("ind must be scalar")
        if (ind <= 0) 
            stop("need positive ind")
        lgd = length(discrim(gd))
        if (ind > length(discrim(gd))) 
            stop(paste("ind request [", ind, "] exceeds number of ",
                 "groups [", lgd, "]"))
        cs = cumsum(as.ram(discrim(gd))) # could be big
        ini = c(1, cs[-lgd]+1)[ind]
        end = cs[ind]
        numdat(gd)[ini:end,]
    }
)
