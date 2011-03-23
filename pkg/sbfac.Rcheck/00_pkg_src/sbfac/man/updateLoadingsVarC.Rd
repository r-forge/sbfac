\name{updateLoadingsVarC}
\alias{updateLoadingsVarC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
insert here the title
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
updateLoadingsVarC(tauinv_, A_, taua_, taub_)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tauinv_}{
%%     ~~Describe \code{tauinv_} here~~
}
  \item{A_}{
%%     ~~Describe \code{A_} here~~
}
  \item{taua_}{
%%     ~~Describe \code{taua_} here~~
}
  \item{taub_}{
%%     ~~Describe \code{taub_} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tauinv_, A_, taua_, taub_) 
.Call("updateLoadingsVarC", tauinv_, A_, taua_, taub_, PACKAGE = "sbfac")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
