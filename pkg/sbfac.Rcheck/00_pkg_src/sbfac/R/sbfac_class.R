sbfac <- function(D, ...) UseMethod("sbfac")

sbfac.default <- function(D, num.factor=1, load.var.type="byrow", 
                          priors = list(taua=0.5, taub=0.5, rhoa=0.5, rhob=0.5), 
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
				obslabel=rownames(D), obslabel=colnames(D),
				nsim=0, nburn=0, thin=1, 
        
        loadings=matrix(rnorm(K*P), nrow=P),
        post.loadings=matrix(rep(0,K*P), ncol=K), 
        post.loadings.var=matrix(rep(0,K*P), ncol=K),
        
        scores=matrix(rnorm(K*N), nrow=K),
        post.scores=matrix(rep(0,K*N), nrow=K),
        post.scores.var=matrix(rep(0,K*N), nrow=K),
        
        tauinv = rep(0, P), rho = rep(0, K),
        priors=priors)
  
	class(out) = "sbfac"
  return(out)
}

print.sbfac <- function(model, ...) {
	cat("SBFAC model object\n")
	cat("Parameter estimates based on",(model$nsim-model$nburn)/model$thin,"MCMC iterations\n")
	cat(model$nsim,"iterations,",model$nburn,"burn-in, keep every",model$thin )
  cat("\nUse summary() for detailed output.")
	}


