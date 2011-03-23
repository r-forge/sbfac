\name{plotLoadingsHeat}
\title{Display a heatmap of factor loadings (note that this function requires ggplot2 package to be...}
\usage{plotLoadingsHeat(model, xlabel=NA, ylabel=NA, 
color=c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", 
"#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"), sorting=NA)}
\description{Display a heatmap of factor loadings (note that this function requires ggplot2 package to be
installed)}
\alias{plotLoadingsHeat}
\arguments{\item{model}{an sbfac model object}
\item{xlabel}{x-axis label}
\item{ylabel}{y-axis label}
\item{color}{character vector of colors}
\item{sorting}{a permutaion of 1:P (where P is the number of variables) providing sort order
for the rows of the loadings matrix}}
