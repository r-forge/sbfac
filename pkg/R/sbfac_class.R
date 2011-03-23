#sbfac <- function(D, ...) UseMethod("sbfac")

#' Initialize an SBFAC model
#'
#' This function accepts a data matrix \code{D} and specified options,
#' returning an S3 object of class sbfac.
#'
#' @param D a numeric matrix where rows are variables and columns are observations
#' @param num.factor number of factors
#' @param load.var.type Not yet implemented
#' @param priors A list of prior parameters. taua and taub are the shape and rate parameters
#' for the gamma prior on the loadings precision tau (i.e. E[tau] = taua/taub) 
#' and rhoa, rhob are the parameters of the beta prior on rho.
#' @param obslabel a character vector of labels for each observation
#' @param varlabel a character vector of labels for each variable
#' @return An S3 object of class \code{sbfac}, which is a list. The elements
#' ldata, loadings, scores, tauinv and rho are initialized values of model parameters. By default
#' these are randomly initialized; you can overwrite them if desired.
#' @export


sbfacModel <- function(D, num.factor=1, load.var.type="byrow", 
                          priors = list(taua=0.5, taub=0.5, rhoa=1.0, rhob=1.0), 
                          obslabel=colnames(D), varlabel=rownames(D), 
                          ...) {
	D = as.matrix(D)
	D[is.infinite(D)] = NA
	p = dim(D)[1]
	n = dim(D)[2]
	Ra<-NULL
	argsorts <- NULL
	for (i in 1:p) {
		Ra<-rbind(Ra, match(D[i,],sort(unique(D[i,]))))
	}
	Ra[is.na(Ra)] = 0
	Ra = Ra - 1 # zero indexing for Rcpp
	argsorts = t(apply(Ra, 1, order))-1

	ecdf.ct <- apply(D,1,rank,ties.method="max",na.last="keep")
	num.not.na = apply(!is.na(D),1,sum)
	U <- t(ecdf.ct)/(num.not.na+1)
	Z = qnorm(U)

	maxes = matrix(rep(0, p*n), nrow=p)
	L = 0
	for (i in 1:p) {
		m = c(sort(unique(Z[i,]))[-1], 7.4)
		maxes[i,1:length(m)] = m
		L = max(L, length(m))
		}
	
	maxes = maxes[,1:L]

  Z[is.na(Z)] = rnorm(length(Z[is.na(Z)]))

	out = list(data=D, ldata=Z, ranks=Ra, maxes=maxes, 	
				argsorts=argsorts, P = dim(D)[1], N = dim(D)[2], K=num.factor,
				obslabel=obslabel, varlabel=varlabel,
				nsim=0, nburn=0, thin=1, 
        
        loadings = matrix(rnorm(K*P), nrow=P),
        post.loadings.mean = matrix(rep(0,K*P), ncol=K), 
        post.loadings.var=matrix(rep(0,K*P), ncol=K),
        
        scores=matrix(rnorm(K*N), nrow=K),
        post.scores.mean=matrix(rep(0,K*N), nrow=K),
        post.scores.var=matrix(rep(0,K*N), nrow=K),
        
        tauinv = rgamma(P, priors$taua)/priors$taub, 
        rho = rbeta(K, priors$rhoa, priors$rhob),
        priors=priors)
  
	class(out) = "sbfac"
  return(out)
}


#' @param x An sbfac object
#' @param ... Ignored
#' @method print sbfac

print.sbfac <- function(x, ...) {
	cat("SBFAC model object\n")
	cat("Parameter estimates based on",x$nsim/x$thin,"MCMC iterations\n")
	cat(x$nsim,"iterations, after",x$nburn,"burn-in, keeping every",x$thin,"samples" )
  cat("\nUse summary() for detailed output.\n")
	}

#' @param x An sbfac object
#' @param ... Ignored
#' @return A list with elements loadings and scores containing MCMC means
#' @method mean sbfac
mean.sbfac <- function(x, ...) {
  return(list(loadings=x$post.loadings.mean, scores=x$post.scores.mean))
}

#' @param x An sbfac object
#' @param ... Ignored
#' @method var sbfac
#' @return A list with elements loadings and scores containing MCMC variances
var.sbfac <- function(x, ...) {
  return(list(loadings=x$post.loadings.var, scores=x$post.scores.var))
}

