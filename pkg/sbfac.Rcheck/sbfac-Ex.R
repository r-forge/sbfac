pkgname <- "sbfac"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('sbfac')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("plotLoadingsHeat")
### * plotLoadingsHeat

flush(stderr()); flush(stdout())

### Name: plotLoadingsHeat
### Title: insert here the title
### Aliases: plotLoadingsHeat
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Z_, A_, F_, tauinv_, rho_) 
.Call("updateSparseLoadingsJ", Z_, A_, F_, tauinv_, rho_, PACKAGE = "sbfac")



cleanEx()
nameEx("updateLoadingsVarC")
### * updateLoadingsVarC

flush(stderr()); flush(stdout())

### Name: updateLoadingsVarC
### Title: insert here the title
### Aliases: updateLoadingsVarC
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tauinv_, A_, taua_, taub_) 
.Call("updateLoadingsVarC", tauinv_, A_, taua_, taub_, PACKAGE = "sbfac")



cleanEx()
nameEx("updateLoadingsVarJC")
### * updateLoadingsVarJC

flush(stderr()); flush(stdout())

### Name: pdateLoadingsVarJC
### Title: insert here the title
### Aliases: pdateLoadingsVarJC
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tauinv_, A_, taua_, taub_) 
.Call("pdateLoadingsVarJC", tauinv_, A_, taua_, taub_, PACKAGE = "sbfac")



cleanEx()
nameEx("updateRho")
### * updateRho

flush(stderr()); flush(stdout())

### Name: updateRho
### Title: insert here the title
### Aliases: updateRho
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rho_, A_, rhoa_, rhob_) 
.Call("updateRho", rho_, A_, rhoa_, rhob_, PACKAGE = "sbfac")



cleanEx()
nameEx("updateScoresC")
### * updateScoresC

flush(stderr()); flush(stdout())

### Name: updateScoresC
### Title: insert here the title
### Aliases: updateScoresC
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Z_, A_, F_) 
.Call("updateScoresC", Z_, A_, F_, PACKAGE = "sbfac")



cleanEx()
nameEx("updateSparseLoadings")
### * updateSparseLoadings

flush(stderr()); flush(stdout())

### Name: updateSparseLoadings
### Title: insert here the title
### Aliases: updateSparseLoadings
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Z_, A_, F_, tauinv_, rho_) 
.Call("updateSparseLoadings", Z_, A_, F_, tauinv_, rho_, PACKAGE = "sbfac")



cleanEx()
nameEx("updateSparseLoadingsJ")
### * updateSparseLoadingsJ

flush(stderr()); flush(stdout())

### Name: updateSparseLoadingsJ
### Title: insert here the title
### Aliases: updateSparseLoadingsJ
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Z_, A_, F_, tauinv_, rho_) 
.Call("updateSparseLoadingsJ", Z_, A_, F_, tauinv_, rho_, PACKAGE = "sbfac")



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
