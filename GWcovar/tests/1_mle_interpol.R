# script to test wether the calculation of the MLE using 'mle.nomean()'
# from the package 'spam' works.

set.seed(43)

require(spam)
require(foreach)
require(mvtnorm)
require(GWcovar)



nbr.col = 5 
bet = 0.5

x<-seq(0,1,len = nbr.col )
loc<-expand.grid(x,x)



###########################
# spam matrices
###########################
dist.mat <- nearest.dist(loc, delta=bet*2, upper=NULL)

result1.2 <- foreach( kap=seq(0,5,by=0.25), .combine='cbind' ) %do% {

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

			covar = cov.wend.interpol( 
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
							  cov.wend.interpol,
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
		sum( result1.2 == -1 ),
		length( result1.2 ),
		"cases mle.nomean() produced an error"
) 
sprintf("[spam matrices 1] In %d of %d %s",
		sum( result1.2 == 0 ),
		length( result1.2 ),
		"cases mle.nomean() did converge"
)
###########################
# dense matrices
###########################

dist.mat <- as.matrix(dist(loc, upper=NULL))

result1.3 <- foreach( kap=seq(0,5,by=0.25), .combine='cbind' ) %do% {

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

			covar = cov.wend.interpol( 
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
							  cov.wend.interpol,
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
		sum( result1.3 == -1 ),
		length( result1.3 ),
		"cases mle.nomean() produced an error"
)
sprintf("[dense matrices 1] In %d of %d %s",
		sum( result1.3 == 0 ),
		length( result1.3 ),
		"cases mle.nomean() did converge"
)


tol2 <- 0.02*length(result1.2)
tol3 <- 0.02*length(result1.3)
if ( sum( result1.2 != 0 ) > tol2 || sum( result1.3 != 0 ) > tol3 ) { 
	string1 = "cases mle.nomean() exited with an error"
	string2 = "cases mle.nomean() did not converge (including errors)"
	error_message = sprintf( 
							rep("\n%s In %d of %d %s",4) ,
							"[spam matrices 0] ",
							sum( result1.2 == -1 ),
							length( result1.2),
							string1,
							"[spam matrices 1] ",
							sum( result1.2 != 0 ),
							length( result1.2 ),
							string2,
							"[dense matrices 0] ",
							sum( result1.3 == -1 ),
							length( result1.3),
							string1,
							"[dense matrices 1] ",
							sum( result1.3 != 0 ),
							length( result1.3 ),
							string2
	)
	stop( error_message )
}
