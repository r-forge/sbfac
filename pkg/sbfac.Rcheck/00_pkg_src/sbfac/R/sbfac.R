 updateScoresC <- function (Z_, A_, F_) 
.Call("updateScoresC", Z_, A_, F_, PACKAGE = "sbfac")

 updateLoadingsVarC <- function (tauinv_, A_, taua_, taub_) 
.Call("updateLoadingsVarC", tauinv_, A_, taua_, taub_, PACKAGE = "sbfac")

 updateLoadingsVarJC <- function (tauinv_, A_, taua_, taub_) 
.Call("updateLoadingsVarJC", tauinv_, A_, taua_, taub_, PACKAGE = "sbfac")

 updateRho <- function (rho_, A_, rhoa_, rhob_) 
.Call("updateRho", rho_, A_, rhoa_, rhob_, PACKAGE = "sbfac")

 updateSparseLoadings <- function (Z_, A_, F_, tauinv_, rho_) 
.Call("updateSparseLoadings", Z_, A_, F_, tauinv_, rho_, PACKAGE = "sbfac")

 updateSparseLoadingsJ <- function (Z_, A_, F_, tauinv_, rho_) 
.Call("updateSparseLoadingsJ", Z_, A_, F_, tauinv_, rho_, PACKAGE = "sbfac")

MCMCstep <- function (Z_, A_, F_, tauinv_, rho_)
.Call("MCMCstep", Z_, A_, F_, tauinv_, rho_,  PACKAGE = "sbfac")