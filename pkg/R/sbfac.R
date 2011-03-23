#' @nord
.updateScoresC <- function (Z_, A_, F_) 
.Call("updateScoresC", Z_, A_, F_, PACKAGE = "sbfac")

#' @nord
 .updateLoadingsVarC <- function (tauinv_, A_, taua_, taub_) 
.Call("updateLoadingsVarC", tauinv_, A_, taua_, taub_, PACKAGE = "sbfac")

#' @nord
 .updateLoadingsVarJC <- function (tauinv_, A_, taua_, taub_) 
.Call("updateLoadingsVarJC", tauinv_, A_, taua_, taub_, PACKAGE = "sbfac")

#' @nord
 .updateRho <- function (rho_, A_, rhoa_, rhob_) 
.Call("updateRho", rho_, A_, rhoa_, rhob_, PACKAGE = "sbfac")

#' @nord
 .updateSparseLoadings <- function (Z_, A_, F_, tauinv_, rho_) 
.Call("updateSparseLoadings", Z_, A_, F_, tauinv_, rho_, PACKAGE = "sbfac")

#' @nord
 .updateSparseLoadingsJ <- function (Z_, A_, F_, tauinv_, rho_) 
.Call("updateSparseLoadingsJ", Z_, A_, F_, tauinv_, rho_, PACKAGE = "sbfac")

#' @nord
#SEXP MCMCstep( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ , SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP priors_, SEXP nsim_, SEXP nburn_, SEXP thin_)
.MCMCstep <- function (Z_, A_, F_, tauinv_, rho_, Ra_, maxes_, argsorts_, priors_, nsim_, nburn_, thin_, printstatus_, keep.scores, keep.loadings)
.Call("MCMCstep", Z_, A_, F_, tauinv_, rho_, Ra_, maxes_, argsorts_, priors_, nsim_, nburn_, thin_, printstatus_, keep.scores, keep.loadings, PACKAGE = "sbfac")

#' Perform MCMC model fitting for an SBFAC model
#'
#' This function performs a specified number of MCMC iterations and
#' returns an sbfac object containing summary statistics from the MCMC samples
#' as well as the actual samples if keep.scores or keep.loadings are TRUE.
#' Default behavior is to save only the loadings. It is recommended to examine
#' traces and marginal posterior density estimates for the loadings as these can be highly skewed
#' and/or multimodal so that the mean/variance are poor summaries. The scores are generally
#' more 'normal'-looking. Take care with these settings as the samples can be very high dimensional.
#'
#' @param model an object of type sbfac, as returned by sbfac(data)
#' @param nsim number of iterations past burn-in
#' @param nburn number of initial (burn-in) iterations to discard
#' @param thin keep every thin'th MCMC sample (i.e. save nsim/thin samples)
#' @param print.status how often to print status messages to console
#' @param keep.scores save samples of factor scores
#' @param keep.loadings save samples of factor loadings
#' @return The S3 \code{sbfac} object \code{model}, now with posterior samples/summaries.
#' @export

doMCMC <- function(model, nsim, nburn, thin=1, print.status=200, keep.scores=FALSE, keep.loadings=TRUE) {
	model$nsim = nsim
	model$nburn = nburn
	model$thin = thin
	sim = .MCMCstep(model$ldata, model$loadings, model$scores, model$tauinv, 
		 model$rho, model$ranks, model$maxes, model$argsorts, 
		 model$priors, nsim, nburn, thin, print.status, keep.scores, keep.loadings)
  
  K = model$K; P = model$P; N = model$N
	
  if (keep.scores)   { dim(sim$Fp)=c(K,N,nsim/thin) }
	if (keep.loadings) { dim(sim$Ap)=c(P,K,nsim/thin) }
	
	dim(sim$Fp.mean) = c(K,N)
	dim(sim$Fp.var)  = c(K,N)
	dim(sim$Ap.mean) = c(P,K)
	dim(sim$Ap.var)  = c(P,K)
	
	model$post.loadings.mean = sim$Ap.mean
	model$post.loadings.var = sim$Ap.var
	model$post.loadings = sim$Ap
	
	model$post.scores.mean = sim$Fp.mean
	model$post.scores.var = sim$Fp.var
	model$post.scores = sim$Fp
	
	return(model)
}