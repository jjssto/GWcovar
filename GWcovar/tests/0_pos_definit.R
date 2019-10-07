# Tests if the calculated covariance matrices are positiv definite

set.seed(43) 

require('spam')
require('foreach')
require('GWcovar')

###########################
# spam matrices
###########################
nbr.col = 5
bet = 0.5

x<-seq(0,1,len = nbr.col )
loc<-expand.grid(x,x)
dist.mat <- nearest.dist(loc, delta=bet, upper=NULL)

result0.0 <- foreach( kap=seq(0,5,by=0.1), .combine='cbind' ) %do% {

	foreach( mu=seq(4.5,9.5,by=0.1 ), .combine='c' ) %do% {

		covar = cov.wend( 
						 dist.mat,
						 c(bet,mu,kap)
		)
		min(eigen(covar)$values)
	}
}

sprintf(
		"[spam matrices] %d of %d matrices are positive definite", 
		sum( result0.0>0 ) ,
		length(result0.0)
)



###########################
# dense matrices
###########################
nbr.col = 5
bet = 0.5

x <- seq(0,1,len = nbr.col )
loc <- expand.grid(x,x)
dist.mat <- as.matrix(dist(loc, upper=NULL))

result0.1 <- foreach( kap=seq(0,5,by=0.1), .combine='cbind' ) %do% {

	foreach( mu=seq(4.5,9.5,by=0.1 ), .combine='c' ) %do% {

		covar = cov.wend( 
						 dist.mat,
						 c(bet,mu,kap)
		)
		min(eigen(covar)$values)
	}
}


sprintf(
		"[dense matrices] %d of %d matrices are positive definite", 
		sum( result0.1>0 ) ,
		length(result0.1)
)

if ( sum( result0.0 <= 0 ) > 0 || sum( result0.1 <= 0 ) > 0 ) { 
	string = "matrices are not positiv definit"
	error_message = sprintf( 
							"\n%s %d of %d %s\n%s %d of %d %s\n" ,
							"[spam matrices] ",
							sum( result0.0 <= 0 ),
							length( result0.0),
							string,
							"[dense matrices] ",
							sum( result0.1 <= 0 ),
							length( result0.1 ),
							string
	)
	stop( error_message )
}

