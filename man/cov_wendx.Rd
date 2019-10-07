\name{cov_wendx}
\alias{cov_wendx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculatete the GW covariance matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cov_wendx(h, theta, abstol = 1e-05, reltol = 0.001, integral_method = 0, qag_intervals = 100, qag_key = 4, n_interpol = 100, fixed_n = 100, eps = getOption("spam.eps"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
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
function (h, theta, abstol = 1e-05, reltol = 0.001, eps = getOption("spam.eps")) 
{
    if ((abstol <= 0) || (abstol <= 0) || (integral_method < 
        0) || (qag_intervals <= 0) || (qag_key < 1) || (qag_key > 
        6) || (n_interpol <= 0)) {
        stop("Invalid parameteres")
    }
    if (length(theta) == 1) {
        theta[2] <- 5.5
        theta[3] <- 1
        theta[4] <- 1
        theta[5] <- 0
    }
    else if (length(theta) == 2) {
        theta[3] <- theta[2]
        theta[2] <- theta[2] + 4.5
        theta[4] <- 1
        theta[5] <- 0
    }
    else if (length(theta) == 3) {
        theta[4] <- 1
        theta[5] <- 0
    }
    else if (length(theta) == 4) {
        theta[5] <- 0
    }
    if (theta[3] >= theta[2]) {
        stop("Kappa bigger than mu.")
    }
    if (is.spam(h)) {
        tryCatch({
            h@entries <- .Call("covar_vector", h@entries, length(h@entries), 
                theta[2], theta[3], theta[4], theta[1], theta[5], 
                abstol, reltol, eps, as.integer(integral_method), 
                as.integer(qag_intervals), as.integer(qag_key), 
                as.integer(n_interpol))
        }, error = function(e) {
            stop("Error")
        })
        return(h)
    }
    else {
        if (integral_method != 1) {
            ret <- .Call("covar_mat", h, theta[2], theta[3], 
                theta[4], theta[1], theta[5], abstol, reltol, 
                as.integer(integral_method), as.integer(n_interpol))
            if (is.null(h)) {
                stop("Invalid matrix")
            }
            else {
                return(ret)
            }
        }
        else {
            stop("Invalid integration method")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
