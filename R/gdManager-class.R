### =========================================================================
### gdManager objects
### -------------------------------------------------------------------------
###

setOldClass("ff_vector")
setOldClass("ff_matrix")
setOldClass("formula")

setClass("gdManager", 
    representation(
        discrim="ff_vector",
        numdat="ff_matrix", 
        formula="formula"
    )
)

### getters / setters
discrim <- function(x) x@discrim
`discrim<-` <- function(x, value) {
    x@discrim <- value
    x
}
numdat <- function(x) x@numdat
`numdat<-` <- function(x, value) {
    x@numdat <- value
    x
}
setMethod("formula", "gdManager", function(x, ...) x@formula) 

### show 
setMethod("show", "gdManager", 
    function(object) {
        cat("gdManager instance\n")
        cat("number of clusters = ", length(object@discrim), "\n")
        nr = nrow(numdat(object))
        nc = ncol(numdat(object))
        cat("size of data matrix: ", nr, "x", nc, "\n")
        cat("excerpt:\n")
        print(numdat(object)[1:min(c(3,nr)),1:min(c(5,nc))])
    }
)

###  groupedData2ff 
groupedData2ff <- function(gd, gcol, prefix="gdm", overwrite=TRUE) 
{
    discrim <- ff(runLength(Rle(gd[,gcol])), 
                  filename=paste(prefix, "disc.ff", sep="."))
    numdat <- ff(data.matrix(gd[,-gcol]), dim=c(nrow(gd), ncol(gd)-1), 
                 dimnames=list(rownames(gd), colnames(gd)[-gcol]),
                 filename=paste(prefix, "dat.ff", sep="."), 
                 overwrite=overwrite)
    new("gdManager", discrim=discrim, numdat=numdat, formula=formula(gd))
}
 
### getGroup 
setGeneric("getGroup", function(x, index, ...) standardGeneric("getGroup")) 
setMethod("getGroup", c("gdManager", "numeric"), 
    function(x, index, ...) {
        if (length(index) > 1) 
            stop("'index' must be scalar")
        if (index <= 0) 
            stop("'index' must be positive")
        ngroup <- length(discrim(x))
        if (index > length(discrim(x))) 
            stop(paste("index request [", index, "] exceeds number of ",
                 "groups [", ngroup, "]"))
        cs <- cumsum(as.ram(discrim(x))) # could be big
        ini <- c(1, cs[-ngroup]+1)[index]
        end <- cs[index]
        numdat(x)[ini:end,]
    }
)

### For unit tests only
setMethod("getGroup", c("ANY", "numeric"), 
    function(x, index, ...) {
        rle <- Rle(getGroups(x))
        groupid <- as.character(getGroupsFormula(x)[[2]])
        cnames <- colnames(x)[!colnames(x) %in% groupid]
        x[as.vector(rle) == runValue(rle)[index], cnames]
    }
)
