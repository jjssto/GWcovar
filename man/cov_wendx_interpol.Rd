\name{cov_wendx_interpol}
\alias{cov_wendx_interpol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculatete the GW covariance matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cov_wendx_interpol(h, theta, kap = 1, abstol = 1e-05, reltol = 0.001, n_interpol = 300, eps = getOption("spam.eps"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{kap}{
%%     ~~Describe \code{kap} here~~
}
  \item{abstol}{
%%     ~~Describe \code{abstol} here~~
}
  \item{reltol}{
%%     ~~Describe \code{reltol} here~~
}
  \item{integral_method}{
%%     ~~Describe \code{integral_method} here~~
}
  \item{qag_intervals}{
%%     ~~Describe \code{qag_intervals} here~~
}
  \item{qag_key}{
%%     ~~Describe \code{qag_key} here~~
}
  \item{n_interpol}{
%%     ~~Describe \code{n_interpol} here~~
}
  \item{fixed_n}{
%%     ~~Describe \code{fixed_n} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
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
function (h, theta, kap = 1, abstol = 1e-05, reltol = 0.001, 
    integral_method = 0, qag_intervals = 100, qag_key = 4, n_interpol = 100, 
    fixed_n = 100, eps = getOption("spam.eps")) 
{
    cov_wendx(h, c(theta[1:2], kap, theta[3:length(theta)]), 
        abstol, reltol, integral_method, qag_intervals, qag_key, 
        n_interpol, fixed_n, eps)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
