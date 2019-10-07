# script to test wether the calculation of the MLE using 'mle.nomean()'
# from the package 'spam' works.

require(spam)
require(foreach)
require(mvtnorm)
require(GWcovar)



nbr.col = 5 
bet = 0.5

x<-seq(0,1,len = nbr.col )
loc<-expand.grid(x,x)



########################################################################
# TEST 1: DOES 'mle.nomean()' produce errors?
########################################################################

###########################
# spam matrices
###########################
dist.mat <- nearest.dist(loc, delta=bet*2, upper=NULL)

result1.0 <- foreach( kap=seq(0,5,by=0.25), .combine='cbind' ) %do% {

	foreach( mu=seq(4.5,7.5,by=0.5 ), .combine='c' ) %do% {

		tryCatch({
			pars <- c(	bet,mu,kap,1,1)
			low <-	c(	1/nbr.col ,
						4.5 ,
						kap/2.0 ,
						0.1,
						0)
			up <-	c( 2*bet,
					   15,
					   15,
					   15,
					   15)

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

result1.1 <- foreach( kap=seq(0,5,by=0.25), .combine='cbind' ) %do% {

	foreach( mu=seq(4.5,7.5,by=0.5 ), .combine='c' ) %do% {

		tryCatch({
			pars <- c(	bet,mu,kap,1,1)
			low <-	c(	1/nbr.col ,
						4.5 ,
						kap/2.0 ,
						0.1,
						0)
			up <-	c( 2*bet,
					   15,
					   15,
					   15,
					   15)

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
		sum( result1.1 == -1 ),
		length( result1.1 ),
		"cases mle.nomean() produced an error"
)
sprintf("[dense matrices 1] In %d of %d %s",
		sum( result1.1 == 0 ),
		length( result1.1 ),
		"cases mle.nomean() did converge"
)


if ( sum( result1.0 != 0 ) > 0 || sum( result1.1 != 0 ) > 0 ) { 
	string1 = "cases mle.nomean() exited with an error"
	string2 = "cases mle.nomean() did not converge (including errors)"
	error_message = sprintf( 
							rep("\nIn %s %d of %d %s",4) ,
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
