\name{sbfacModel}
\alias{sbfacModel}
\title{Initialize an SBFAC model...}
\usage{sbfacModel(D, num.factor=1, load.var.type="byrow", priors=list(taua = 0.5, taub = 0.5, rhoa =
    1, rhob = 1), obslabel=colnames(D), varlabel=rownames(D), ...)}
\description{Initialize an SBFAC model}
\details{This function accepts a data matrix \code{D} and specified options,
returning an S3 object of class \code{sbfac} 
as well as the actual samples if keep.scores or keep.loadings are TRUE.
Default behavior is to save traces of the loadings but not the scores. 
The loadings are more prone to mixing problems.}
\value{An S3 object of class \code{sbfac}, which is a list. The elements
ldata, loadings, scores, tauinv and rho are initialized values of model parameters. By default
these are randomly initialized; you can overwrite them if desired.}
\arguments{\item{D}{a numeric matrix where rows are variables and columns are observations}
\item{num.factor}{number of factors}
\item{load.var.type}{Not yet implemented}
\item{priors}{A list of prior parameters. taua and taub are the shape and rate parameters
for the gamma prior on tau (i.e. E[tau] = taua/taub) and rhoa, rhob are the parameters of the 
beta prior on rho.}
\item{obslabel}{a character vector of labels for each observation}
\item{varlabel}{a character vector of labels for each variable}}
