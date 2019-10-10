# script to test wether the calculation of the MLE using 'mle.nomean()'
# from the package 'spam' works.

set.seed(42)

require(spam)
require(mvtnorm)
require(GWcovar)



nbr.col = 5 
bet = 0.5

x<-seq(0,1,len = nbr.col )
loc<-expand.grid(x,x)

kappas  <- seq(0,5,by=0.25)
mus <- seq(4.5,7.5,by=0.5 )
n_k <- length(kappas)
n_m <- length(mus)

result1.0 <- matrix(NA, n_k, n_m)
result1.1 <- matrix(NA, n_k, n_m)




###########################
# spam matrices
###########################
dist.mat <- nearest.dist(loc, delta=bet*2, upper=NULL)

for( i in 1:n_k ) {
    for( j in 1:n_m ) {

	result1.0[i,j] <- tryCatch({
	    pars <- c(
		      bet,
		      mus[j],
		      kappas[i],
		      1,
		      1
	    )
	    low <- c(
		     1/nbr.col ,
		     4.5 ,
		     kappas[i]/2.0 ,
		     0.1,
		     0
	    )
	    up <- c( 
		    2*bet,
		    15,
		    15,
		    15,
		    15
	    )	   
	    covar = cov.wend( 
			     dist.mat,
			     pars
	    )
	    y <- rmvnorm( 
			 n = 1,
			 sigma = covar
	    )
	    ret <- mle.nomean(
			      t(y),
			      dist.mat,
			      cov.wend,
			      pars, 
			      low,
			      up
	    )	
	    ret$convergence
	}, error = function(e) {
	    return( -1 ) 
	})
    }
}

sprintf("[spam matrices 0] In %d of %d %s",
	sum( result1.0 == -1 ),
	length( result1.0 ),
	"cases mle.nomean() produced an error"
) 

sprintf("[spam matrices 1] In %d of %d %s",
	sum( result1.0 == 0 ),
	length( result1.0 ),
	"cases mle.nomean() did converge"
)


###########################
# dense matrices
###########################

dist.mat <- as.matrix(dist(loc, upper=NULL))

for( i in 1:n_k ) {
    for( j in 1:n_m ) {

	result1.1[i,j] <- tryCatch({
	    pars <- c(
		      bet,
		      mus[j],
		      kappas[i],
		      1,
		      1
	    )
	    low <- c(
		     1/nbr.col ,
		     4.5 ,
		     kappas[i]/2.0 ,
		     0.1,
		     0
	    )
	    up <- c( 
		    2*bet,
		    15,
		    15,
		    15,
		    15
	    )

	    covar = cov.wend( 
			     dist.mat,
			     pars
	    )
	    y <- rmvnorm( 
			 n = 1,
			 sigma = covar
	    )
	    ret <- mle.nomean(
			      t(y),
			      dist.mat,
			      cov.wend,
			      pars, 
			      low,
			      up
	    )	
	    ret$convergence
	}, error = function(e) {
	    return( -1 ) 
	})
    }
}




sprintf("[dense matrices 0] In %d of %d %s",
	sum( result1.1 == -1 ),
	length( result1.1 ),
	"cases mle.nomean() produced an error"
)
sprintf("[dense matrices 1] In %d of %d %s",
	sum( result1.1 == 0 ),
	length( result1.1 ),
	"cases mle.nomean() did converge"
)



#######################
# Evaluation
#######################


# percentage of non-convergence of mle.nomean() that is accepted
tolerance <- 0.01

tol0 <- tolerance*length(result1.0)
tol1 <- tolerance*length(result1.1)


if ( sum( result1.0 != 0 ) > tol0 || sum( result1.1 != 0 ) > tol1 ) { 
    string1 = "cases mle.nomean() exited with an error"
    string2 = "cases mle.nomean() did not converge (including errors)"
    error_message = sprintf( 
			    paste(rep("\n%s In %d of %d %s",4)) ,
			    "[spam matrices 0] ",
			    sum( result1.0 == -1 ),
			    length( result1.0),
			    string1,
			    "[spam matrices 1] ",
			    sum( result1.0 != 0 ),
			    length( result1.0 ),
			    string2,
			    "[dense matrices 0] ",
			    sum( result1.1 == -1 ),
			    length( result1.1),
			    string1,
			    "[dense matrices 1] ",
			    sum( result1.1 != 0 ),
			    length( result1.1 ),
			    string2
    )
    stop( error_message )
}
