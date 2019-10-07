## Copyright (C) 2019 Josef Stocker <josef@josefstocker.ch>
## 
## This program is free software: you can redistribute it and/or modify
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

dyn.load('src/covar.so')


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
