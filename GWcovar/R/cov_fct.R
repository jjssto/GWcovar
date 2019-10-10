## This file is part of the R-package 'GWcovar'
##
## Copyright (C) 2019 Josef Stocker <josef@josefstocker.ch>
## 
## 'GWcovar' is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <https://www.gnu.org/licenses/>.  
##
########################################################################

#' This package provides functions to calculate the Generalized Wendland
#' covariance matrix
#'
#' The package rovides the two functions
#' \code{\link{cov.wend}} and \code{\link{cov.wend.interpol}}.
#' 
#' @author Josef Stocker, \email{josef@@josefstocker.ch}
#' @seealso \pkg{spam}
#' @docType package
#' @name GWcovar
#' @references \url{https://arxiv.org/abs/1607.06921}
#'
#' @import spam 
#' @useDynLib covar, .registration = TRUE

dyn.load('src/covar.so')


#' Calculates the Generalized Wendland covariance matrix.
#'
#' The function \code{cov.wend} calculates the Generalized Wendland (GW)
#' covariance matrix based on a distance matrix. 
#'
#' @return If the distance matrix is in standard R format a standard R matrix is
#' returned. If the distance matrix is of class \linkS4class{spam} the returned matrix is
#' also of class \linkS4class{spam}.
#' 
#' @param h distance matrix
#' @param theta parameter vector (only range range needs to be specified):
#'     theta[1]: range
#'     theta[2]: mu - kappa (default: 5)
#'     theta[3]: kappa (default: 1)
#'     theta[4]: sill (default: 1)
#'     theta[5]: nugget (default: 0)
#' @param abstol absolute tolerance used for the calculation of the GW
#' covariance function
#' @param reltol relative tolerance used for the calculation of the GW
#' covariance function
#' @param eps treshhold below which values are considered to be equal to
#' 0
#'
#' @seealso \linkS4class{spam}}
#' @export
#' @examples
#' x <- seq(0,1,len=10) 
#' loc <- expand.grid(x,x) 
#' dist.mat <- spam::nearest.dist(loc,upper=NULL,delta=0.5)
#' cov.wend( dist.mat, c(0.3,6,1.5,1,0))
cov.wend <- function( 
                      h, 
                      theta, 
                      abstol = 1e-5, 
                      reltol = 1e-2, 
                      eps = getOption("spam.eps")) {

    if ( (abstol <= 0) || (abstol <= 0) || (eps < 0) ) {
        stop("Invalid arguments")
    }
    # Calculates GW covariance function. 
    if(length(theta)==1){
		if ( theta[1] <= 0 ) {
        	stop("Invalid arguments")
		}
        theta[2] <- 5
        theta[3] <- 1.0
        theta[4] <- 1.0
        theta[5] <- 0
    } else if (length(theta)==2) {
		if ( (theta[1]<=0) || (theta[2]<=0) ) {
        	stop("Invalid arguments")
		}
        theta[3] <- 1.5
        theta[4] <- 1.0
        theta[5] <- 0
    } else	if(length(theta)==3) {
		if ( (theta[1]<=0) || (theta[2]<=0) || (theta[3]<0) ) {
        	stop("Invalid arguments")
		}
        theta[4] <- 1.0
        theta[5] <- 0
    } else if ( length(theta)==4 ) {
		if ( (theta[1]<=0) || (theta[2]<=0) || (theta[3]<0) || (theta[4]<=0) ) {
        	stop("Invalid arguments")
		}
        theta[5] <- 0
    } else {
		if ( (theta[1]<=0) || (theta[2]<=0) || (theta[3]<0) || 
			(theta[4]<=0) || (theta[5]<0) ) {
        	stop("Invalid arguments")
		}
	}


    if(spam::is.spam(h)) {

        tryCatch({
            h@entries  <- .Call("covar_vector_dir",	
                                h@entries, length(h@entries), theta[2]+theta[3],theta[3],
                                theta[4],theta[1],theta[5], abstol, reltol, eps
            )
        }, error=function(e) {

			stop("An error occured in the calculation of the covariance matrix.")
        })
        return(h)
    } else {
        ret  <- .Call("covar_m_dist",
                      h ,  theta[2]+theta[3],theta[3],
                      theta[4],theta[1],theta[5], abstol, reltol
        )
        if (is.null(h) ) {

			stop("An error occured in the calculation of the covariance matrix.")
        } else {

            return (ret)
        }
    }
}


#' Calculates the Generalized Wendland covariance matrix.
#'
#' The function \code{cov.wend.interpol} calculates the Generalized
#' Wendland (GW) covariance matrix based on a distance matrix. In
#' contrast to \code{\link{cov.wend}} this function uses interpolation in order to
#' increase calculation speed.
#'
#' @return If the distance matrix is in standard R format a standard R matrix is
#' returned. If the distance matrix is of class 'spam' the returned matrix is
#' also of class \linkS4class{spam}.
#'
#' @param h distance matrix
#' @param theta parameter vector (only range range needs to be specified):
#'     theta[1]: range
#'     theta[2]: mu - kappa (default: 5)
#'     theta[3]: kappa (default: 1)
#'     theta[4]: sill (default: 1)
#'     theta[5]: nugget (default: 0)
#' @param abstol absolute tolerance used for the calculation of the GW
#' covariance function
#' @param reltol relative tolerance used for the calculation of the GW
#' covariance function
#' @param n_interpol number of equidistant locations where the GW
#' covariance function is calculated 
#' @param eps treshhold below which values are considered to be equal to
#' 0
#'
#' @seealso \pkg{spam}
#' @export
#' @examples
#' x <- seq(0,1,len=10) 
#' loc <- expand.grid(x,x) 
#' dist.mat <- spam::nearest.dist(loc,upper=NULL,delta=0.5) 
#' cov.wend.interpol( dist.mat, c(0.3,6,1.5,1,0))
cov.wend.interpol <- function( 
                      h, 
                      theta, 
                      abstol = 1e-5, 
                      reltol = 1e-2, 
                      n_interpol = 300,
                      eps = getOption("spam.eps")) {

    if ( (abstol <= 0) || (abstol <= 0) || (eps<0) || (n_interpol<=0) ) {
        stop("Invalid arguments")
    }
    if(length(theta)==1){
		if ( theta[1] <= 0 ) {
        	stop("Invalid arguments")
		}
        theta[2] <- 5
        theta[3] <- 1
        theta[4] <- 1.0
        theta[5] <- 0
    } else if (length(theta)==2) {
		if ( (theta[1]<=0) || (theta[2]<=0) ) {
        	stop("Invalid arguments")
		}
        theta[3] <- 1
        theta[4] <- 1.0
        theta[5] <- 0
    } else	if(length(theta)==3) {
		if ( (theta[1]<=0) || (theta[2]<=0) || (theta[3]<0) ) {
        	stop("Invalid arguments")
		}
        theta[4] <- 1.0
        theta[5] <- 0
    } else if ( length(theta)==4 ) {
		if ( (theta[1]<=0) || (theta[2]<=0) || (theta[3]<0) || (theta[4]<=0) ) {
        	stop("Invalid arguments")
		}
        theta[5] <- 0
    } else {
		if ( (theta[1]<=0) || (theta[2]<=0) || (theta[3]<0) || 
			(theta[4]<=0) || (theta[5]<0) ) {
        	stop("Invalid arguments")
		}
	}
    if(spam::is.spam(h)) {

		tryCatch({
        	h@entries  <- .Call("covar_vector_interpol",	
                            h@entries, length(h@entries), theta[2]+theta[3],theta[3],
                            theta[4],theta[1],theta[5], abstol, reltol, eps,
                            as.integer(n_interpol) )
		}, error=function(e) {

			stop("An error occured in the calculation of the covariance matrix.")
		})
       	return(h)
    } else {
        ret  <- .Call("covar_interpol",
                      h , theta[2]+theta[3],theta[3],
                      theta[4],theta[1],theta[5], 
                      as.integer(n_interpol),
                      abstol, reltol
        )
        if (is.null(h) ) {

			stop("An error occured in the calculation of the covariance matrix.")
        } else {

            return (ret)
        }
    }
}
