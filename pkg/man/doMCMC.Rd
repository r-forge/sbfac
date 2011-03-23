\name{doMCMC}
\alias{doMCMC}
\title{Perform MCMC model fitting for an SBFAC model...}
\usage{doMCMC(model, nsim, nburn, thin=1, print.status=200, keep.scores=FALSE, keep.loadings=TRUE)}
\description{Perform MCMC model fitting for an SBFAC model}
\details{This function performs a specified number of MCMC iterations and
returns an sbfac object containing summary statistics from the MCMC samples
as well as the actual samples if keep.scores or keep.loadings are TRUE.
Default behavior is to save only the loadings. It is recommended to examine
traces and marginal posterior density estimates for the loadings as these can be highly skewed
and/or multimodal so that the mean/variance are poor summaries. The scores are generally
more 'normal'-looking. Take care with these settings as the samples can be very high dimensional.}
\value{The S3 \code{sbfac} object \code{model}, now with posterior samples/summaries.}
\arguments{\item{model}{an object of type sbfac, as returned by sbfac(data)}
\item{nsim}{number of iterations past burn-in}
\item{nburn}{number of initial (burn-in) iterations to discard}
\item{thin}{keep every thin'th MCMC sample (i.e. save nsim/thin samples)}
\item{print.status}{how often to print status messages to console}
\item{keep.scores}{save samples of factor scores}
\item{keep.loadings}{save samples of factor loadings}}
