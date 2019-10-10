# Tests if the calculated covariance matrices are positiv definite

set.seed(42) 

require('spam')
require('GWcovar')

###########################
# spam matrices
###########################
nbr.col <- 5
bet <- 0.5

x<-seq(0,1,len = nbr.col )
loc<-expand.grid(x,x)

kappas  <- seq(0,5,by=0.1)
mus <- seq(4.5,9.5,by=0.1 )
n_k <- length(kappas)
n_m <- length(mus)

result0.0 <- matrix(NA, n_k, n_m)
result0.1 <- matrix(NA, n_k, n_m)


dist.mat <- nearest.dist(loc, delta=bet, upper=NULL)

for ( i in 1:n_k ) {
    for( j in 1:n_m ) {
	covar = cov.wend( 
			 dist.mat,
			 c(bet,mus[j],kappas[i])
	)
	result0.0[i,j] <- min(eigen(covar)$values) 
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
dist.mat <- as.matrix(dist(loc, upper=NULL))

for ( i in 1:n_k ) {
    for( j in 1:n_m ) {
	covar = cov.wend( 
			 dist.mat,
			 c(bet,mus[j],kappas[i])
	)
	result0.1[i,j] <- min(eigen(covar)$values) 
    }
}


sprintf(
	"[dense matrices] %d of %d matrices are positive definite", 
	sum( result0.1>0 ) ,
	length(result0.1)
)



##############################
# Evaluation
#############################

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

