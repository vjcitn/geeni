### R code from vignette source 'rungee.Rnw'

###################################################
### code chunk number 1: getd
###################################################
library(geeni)
library(nlme)
# from help(groupedData)
  Orth.new <-  # create a new copy of the groupedData object
       groupedData( distance ~ age | Subject,
                   data = as.data.frame( Orthodont ),
                   FUN = mean,
                   outer = ~ Sex,
                   labels = list( x = "Age",
                     y = "Distance from pituitary to pterygomaxillary fissure" ),
                   units = list( x = "(yr)", y = "(mm)") )
dim(Orth.new)
library(ff)
library(geeni)
#targdir = system.file("ffdemo", package="geeni")
#mantargdir = system.file("mgrs", package="geeni")
#pref = paste(targdir, "gdm", sep="/")
#fis = dir(pref, full=TRUE)
#if (length(fis)>0) try(sapply(fis, file.remove))
flatOrth = geeni:::gd2flat(gd=Orth.new, gcol=3, prefix="")# prefix=pref)
flatOrth


###################################################
### code chunk number 2: lkg
###################################################
getGrp(flatOrth,1)
getGrp(flatOrth,4)


###################################################
### code chunk number 3: glminf
###################################################
geeni:::getDep
geeni:::getEta
geeni:::getMu
geeni:::getX 
geeni:::getY
geeni:::Di
geeni:::Vinv.i
geeni:::ri


###################################################
### code chunk number 4: doco
###################################################
beta = c(0,0,0)
delb = function(gd, beta, family) {
 DD = geeni:::Di(gd,1,beta,family)
 val = t(DD) %*% geeni:::Vinv.i(gd,1,beta,family) 
 val1 = val %*% DD
 val2 = val %*% geeni:::ri(gd,1,beta,family)
 for  (i in 2:length(gd@discrim)) {
    DD = geeni:::Di(gd, i, beta, family)
    val = t(DD) %*% geeni:::Vinv.i(gd,i,beta,family) 
    val1 = val1 + val %*% DD
    val2 = val2 + val %*% geeni:::ri(gd,i,beta,family)
 }
 solve(val1)%*%val2
}
delb( flatOrth, beta, gaussian )


###################################################
### code chunk number 5: lkg
###################################################
library(nlme)
example(groupedData)


###################################################
### code chunk number 6: lklm
###################################################
lm(distance~age+Sex,data=Orth.new)


###################################################
### code chunk number 7: lkf
###################################################
geeni:::Gcomps


###################################################
### code chunk number 8: accum
###################################################
geeni:::combi


###################################################
### code chunk number 9: howto
###################################################
library(foreach)
library(doParallel)
registerDoParallel(cores=2)  # for mac
comps = foreach(i = 1:27, .combine=geeni:::combi) %dopar% 
    { geeni:::Gcomps(flatOrth,i,c(0,0,0),gaussian) } 
comps
del = solve(comps[[1]])%*%comps[[2]]
beta = beta + del
comps = foreach(i = 1:27, .combine= geeni:::combi) %dopar% 
    { geeni:::Gcomps(flatOrth,i,beta,gaussian) } 
comps


###################################################
### code chunk number 10: dog
###################################################
geeni::pargee
pargee( flatOrth, gaussian, c(0,0,0) )


###################################################
### code chunk number 11: domore (eval = FALSE)
###################################################
## flatOrth@numdat[,"Sex"] = flatOrth@numdat[,"Sex"]-1  # overwrite allowed


###################################################
### code chunk number 12: killff
###################################################
system("rm -rf .dat.ff")
system("rm -rf .disc.ff")


###################################################
### code chunk number 13: lkses
###################################################
sessionInfo()


