/* Copyright (C) 2019 Josef Stocker <josef@josefstocker.ch>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 * */
#include "covar.h"

#include "stdio.h"
#include "stdlib.h"

#include "R.h"
#include "Rinternals.h"
#include "Rmath.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_errno.h"

#include "wendland.h"


SEXP 
covar_m_dist (  
        SEXP DIST ,      /* distance matrix */
        SEXP MU ,        /* param. of GW correlation function */
        SEXP SMOOTHNESS ,/* param. of GW correlation function */
        SEXP SILL ,      /* param. of GW correlation function */
        SEXP RNGE ,      /* param. of GW correlation function */
        SEXP NUGGET,     /* param. of GW correlation function */
        SEXP ABSTOL,     /* absolute tolerance for integration */
        SEXP RELTOL      /* relative tolerance for integration */
        )
/* **************************************************************************** 
 * The function 'SEXP covar_m_dist(...)' returns the Generalized Wendland
 * (GW) covariance matrix based on a distance matrix in standard R matrix
 * format.  The values of the covariance function are calculated using  
 * non-addaptive Gauss-Kronrod integration.
 * ***************************************************************************/
{
    /* create local representatives for the SEXPs */
    int* p_dim = INTEGER(getAttrib( DIST, R_DimSymbol )) ;
    double* p_dist = REAL( DIST ) ;
    double mu = *REAL( MU ) ;
    double smoothness = * REAL( SMOOTHNESS ) ;
    double sill = *REAL( SILL ) ;
    double rnge = *REAL( RNGE ) ;
    double nugget = *REAL( NUGGET ) ;
    double abstol = *REAL( ABSTOL ) ;
    double reltol = *REAL( RELTOL ) ;

    /* allocate return object */
    SEXP RESULT ;
    PROTECT( 
            RESULT = allocMatrix( REALSXP, *p_dim, *(p_dim+1) ) 
           ) ; 

    Wendland_result result ;
    gsl_set_error_handler_off() ;
    if ( *p_dim == *(p_dim+1) ) {
        /* if matrix is square */

        for ( int i=0 ; i < *p_dim ; i++ ) {
            /* for each row */

            REAL(RESULT)[i + i*(*p_dim)] = sill+nugget ; 
            /* diagonal element  of matrix */

            /* For loop calculates the values for the upper triangular matrix
             * and than uses the simetry of the covariance matrix to insert the
             * value into the lower triangular matrix, too.*/
            for ( int j=i+1 ; j< *(p_dim+1) ; j++ ) {

                if( *(p_dist + i + j * (*p_dim) )  == 0 ) {
                    /* if distance = 0 */

                    REAL(RESULT)[i + j*(*p_dim)] = sill + nugget ;
                    /* upper triangular matrix */

                    REAL(RESULT)[j + i*(*p_dim)] = REAL(RESULT)[i + j*(*p_dim)] ;
                    /* lower triangular matrix */

                } else if (*(p_dist + i + j*(*p_dim) ) < rnge ) {
                    /* if distance < rnge */

                    /* upper triangular matrix */
                     wendland( &result,
                            *(p_dist+i+j * (*p_dim) ) / rnge ,
                            mu, smoothness, abstol, reltol ) ;


                    if ( check_wendland_errors( &result ) ) {

                        REAL(RESULT)[i + j*(*p_dim)] = sill * result.result ;

                        REAL(RESULT)[j + i*(*p_dim)] = REAL(RESULT)[i + j*(*p_dim)] ;
                    } else {

                        return R_NilValue ;
                    }



                    /* lower triangular matrix */
                } else {

                    REAL(RESULT)[i + j*(*p_dim)] = 0 ;
                    /* upper triangular matrix */

                    REAL(RESULT)[j + i*(*p_dim)] = 0 ;
                    /* lower triangular matrix */
                }
            }
        }
    } else {
        /* if matrix is not square */

        /* For loop iterates through all matrix entries */
        for ( int i=0 ; i< *p_dim  ; i++ ) {

            for ( int j=0; j< *(p_dim+1) ; j++ ) {

                if( *(p_dist + i + j * (*p_dim) ) == 0 ) {

                    REAL(RESULT)[i + j*(*p_dim)] = sill + nugget ;
                } else if ( *(p_dist+i + j * (*p_dim) ) < rnge ) {

                    wendland( &result ,
                            *(p_dist +i +j * (*p_dim) ) / rnge, 
                            mu, smoothness, abstol, reltol ) ;
                    
                    if ( check_wendland_errors( &result ) ) {

                        REAL(RESULT)[i + j*(*p_dim)] = sill * result.result ;
                    } else {

                        return R_NilValue ;
                    }
                } else {

                    REAL(RESULT)[i + j*(*p_dim)] = 0 ;
                }
            }

        }
    }
    UNPROTECT(1) ; /* RESULT */
    return RESULT ;
}




SEXP 
covar_interpol (
        SEXP DIST ,     /* distance matrix */
        SEXP MU ,       /* param. of GW correlation function */
        SEXP SMOOTHNESS,/* param. of GW correlation function */
        SEXP SILL ,     /* param. of GW correlation function */
        SEXP RNGE ,     /* param. of GW correlation function */
        SEXP NUGGET ,   /* param. of GW correlation function */
        SEXP N ,        /* nbr. of points where fct is
                         * evaluated for interpolation */ 
        SEXP ABSTOL ,   /* rel. tolerance for integration */
        SEXP RELTOL     /* abs. tolerance for integration */ 
        )
/* **************************************************************************** 
 * The function 'SEXP covar_mat_interpol(...)' returns the Generalized Wendland
 * (GW) covariance matrix based on a distance matrix in standard R matrix
 * format.  The values of the covariance function are interpolated using cubic
 * splines. For the calculation of the GW covariance function the
 * non-addaptive Gauss-Kronrod algorithm is used.
 * ***************************************************************************/
{
    /* create local representation for the SEXPs */
    int* p_dim = INTEGER(getAttrib( DIST, R_DimSymbol )) ;
    double* p_dist = REAL( DIST ) ;
    double mu = *REAL( MU ) ;
    double smoothness = * REAL( SMOOTHNESS ) ;
    double sill = *REAL( SILL ) ;
    double rnge = *REAL( RNGE ) ;
    double nugget = *REAL( NUGGET ) ;
    int n = *INTEGER( N ) ;
    double abstol = *REAL( ABSTOL ) ;
    double reltol = *REAL( RELTOL ) ;

    
    /* allocate return object */
    SEXP RESULT ;
    PROTECT( 
            RESULT = allocMatrix( REALSXP, *p_dim, *(p_dim+1) )
           ) ; 

    double *points = malloc( n * sizeof(double) ) ;
    /* allocating array for the interpolation points */

    double *wendl = malloc( n * sizeof(double) ) ;
    /* allocation array for the values of the cov. fct.  in the interpolation
     * points */

    double interval = rnge / ( (double) n - 1.0 ) ;
    /* distance between the interpolation points */

    gsl_set_error_handler_off() ;
    Wendland_result result ;

    /* calculating the function in the interpolation points */
    for ( int i=0 ; i < n ; ++i ) {

        points[i] = i*interval ;
        wendland( &result,  points[i]/rnge, mu, smoothness, abstol, reltol ) ;

        if ( check_wendland_errors( &result ) ) {

            wendl[i]    = sill * result.result ;     
        } else {

            return R_NilValue ;
        }   
    }

    /* initialisation for the interpolation */
    gsl_interp *interpol = gsl_interp_alloc(gsl_interp_cspline , n ) ;
    gsl_interp_init( interpol,  points, wendl, n ) ;
    gsl_interp_accel *acc =  gsl_interp_accel_alloc() ;

    if ( *p_dim == *(p_dim+1) ) {
        /* if the matrix is square */

        /* For loop iterates through the upper triangular matrix. For the lower
         * triangular matrix symetry is used. */
        for ( int i=0 ; i< *p_dim ; i++ ) {

            REAL(RESULT)[i + i*(*p_dim)] = sill+nugget ; 
            /* diagonal of matrix */

            for ( int j=i+1 ; j< *(p_dim+1) ; j++ ) {

                if ( *(p_dist+i+j*(*p_dim)) == 0 ) {
                    /* dist = 0 */ 

                    REAL(RESULT)[i + j*(*p_dim)] = sill + nugget ;
                    /* upper triangular matrix */

                    /* lower triangular matrix */
                    REAL(RESULT)[j + i*(*p_dim)] =
                        REAL(RESULT)[i + j*(*p_dim)] ;
                } else if ( *(p_dist+i+j*(*p_dim)) < rnge ) {
                    /* dist < rnge */

                    /* upper triangular matrix */
                    REAL(RESULT)[i + j*(*p_dim)] = 
                        gsl_interp_eval( interpol, points, wendl, 
                                *(p_dist+i+j*(*p_dim)) , acc ) ;

                    /* lower triangular matrix */
                    REAL(RESULT)[j + i*(*p_dim)] =
                        REAL(RESULT)[i + j*(*p_dim)] ;
                } else {
                    /* dist > rnge */

                    REAL(RESULT)[i + j*(*p_dim)] = 0 ;
                    /* upper triangular matrix */

                    REAL(RESULT)[j + i*(*p_dim)] = 0 ;
                    /* lower triangular matrix */
                }
            }
        }
    } else {
        /* matrix is not square */

        /* For loop iterates through all matrix entries */
        for ( int i=0 ; i< *p_dim ; i++ ) {

            for ( int j=0; j< *(p_dim+1) ; j++ ) {

                if( *(p_dist + i + j * (*p_dim)) == 0 ) {
                    /* dist = 0 */

                    REAL(RESULT)[i + j*(*p_dim)] = sill + nugget ;
                } else if ( *(p_dist+i+j*(*p_dim)) < rnge ) {
                    /* dist < rnge */

                    REAL(RESULT)[i + j*(*p_dim)] = gsl_interp_eval( 
                            interpol, points, wendl, *(p_dist+i+j*(*p_dim)) , acc ) ;
                } else {
                    /* dist > rnge */

                    REAL(RESULT)[i + j*(*p_dim)] = 0 ;
                }
            }

        }
    }
    gsl_interp_accel_free( acc ) ;
    gsl_interp_free( interpol ) ;   
    free(points) ;
    free(wendl) ;
    UNPROTECT(1) ; /* RESULT */
    return RESULT ;
}

SEXP covar_vector_dir (
        SEXP DIST ,         /* R vector containing distances */    
        SEXP LENGTH ,       /* length of 'SEXP DIST' */ 
        SEXP MU ,           /* param. of the GW covariance fct */
        SEXP SMOOTHNESS ,   /* param. of the GW covariance fct */
        SEXP SILL ,         /* param. of the GW covariance fct */
        SEXP RNGE ,         /* param. of the GW covariance fct */
        SEXP NUGGET ,       /* param. of the GW covariance fct */
        SEXP ABSTOL ,       /* abs. tolerance for integration */
        SEXP RELTOL ,       /* rel. tolerance for integration */
        SEXP EPS            /* treshhold below which values are
                             * considered 0 */
       )
/* ****************************************************************************
* The function 'int covar_vector_dir  (...)' calculates the Generalized Wendland
* (GW) covariance matrix by modifying the R vector 'DIST'. The function
* returns 0 if the calculation finished without erros.  Depending on the
* parameters used the integral is calculated with the non-adaptive
* Gauss-Kronrod algorithm using the 'GNU Scientific Library'. 
* **************************************************************************/
{
    /* local representation for the SEXPs */
    double* p_dist = REAL(DIST) ;
    int length = *INTEGER( LENGTH ) ;
    double mu = *REAL( MU ) ;
    double smoothness = *REAL( SMOOTHNESS ) ;
    double sill = *REAL( SILL ) ;
    double rnge = *REAL( RNGE ) ;
    double nugget = *REAL( NUGGET ) ;
    double abstol = *REAL( ABSTOL ) ;
    double reltol = *REAL( RELTOL ) ;
    double eps = *REAL( EPS ) ;

    /* declare and allocate matrix that will be returned */
    SEXP RESULT ;
    PROTECT( 
            RESULT = allocVector( REALSXP, length ) 
           ) ;

    gsl_set_error_handler_off() ;
    Wendland_result result ;

    for( int i = 0 ; i < length ; i++ ) {

        if ( *(p_dist+i) < eps ) {

            REAL(RESULT)[i] = sill + nugget ;
        } else {


            wendland ( 

                    &result,
                    *(p_dist+i)/rnge,  
                    mu, 
                    smoothness, 
                    abstol, 
                    reltol
                    ) ;

            if ( check_wendland_errors( &result ) ) {

                    REAL(RESULT)[i] = sill * result.result ;
            } else {

                return R_NilValue ;
            }
        }         /* if clause */
    } /* for loop */

    UNPROTECT(1) ; /* RESULT */
    return RESULT ;
}

SEXP covar_vector_interpol (
        SEXP DIST ,         /* R vector containing distances */    
        SEXP LENGTH ,       /* length of 'SEXP DIST' */ 
        SEXP MU ,           /* param. of the GW covariance fct */
        SEXP SMOOTHNESS ,   /* param. of the GW covariance fct */
        SEXP SILL ,         /* param. of the GW covariance fct */
        SEXP RNGE ,         /* param. of the GW covariance fct */
        SEXP NUGGET ,       /* param. of the GW covariance fct */
        SEXP ABSTOL ,       /* abs. tolerance for integration */
        SEXP RELTOL ,       /* rel. tolerance for integration */
        SEXP EPS ,          /* treshhold below which values are
                             * considered 0 */
        SEXP NBR_INTERPOL  /* nbr. of interpolation points */ 
        )
/* *****************************************************************************
 * The function 'int covar_vector_dir  (...)' calculates the Generalized Wendland
 * (GW) covariance matrix by modifying the R vector 'DIST'. The function
 * returns 0 if the calculation finished without erros.  Depending on the
 * parameters used the integral is calculated with the non-adaptive
 * Gauss-Kronrod algorithm using the 'GNU Scientific Library'. 
 * **************************************************************************/
{
    /* local representation for the SEXPs */
    double* p_dist = REAL(DIST) ;
    int length = *INTEGER( LENGTH ) ;
    double mu = *REAL( MU ) ;
    double smoothness = *REAL( SMOOTHNESS ) ;
    double sill = *REAL( SILL ) ;
    double rnge = *REAL( RNGE ) ;
    double nugget = *REAL( NUGGET ) ;
    double abstol = *REAL( ABSTOL ) ;
    double reltol = *REAL( RELTOL ) ;
    double eps = *REAL( EPS ) ;
    int n = *INTEGER( NBR_INTERPOL ) ;


    /* declare and allocate matrix that will be returned */
    SEXP RESULT ;
    PROTECT( 
            RESULT = allocVector( REALSXP, length ) 
           ) ;


    /* interpolation */

    double *points = malloc( n * sizeof(double) ) ;
    /* allocating array for the interpolation points */

    double *wendl = malloc( n * sizeof(double) ) ;
    /* allocatin array for the covariance values in the
     * interpolation points */

    double interval = rnge / ( (double) n - 1.0 ) ;
    /* distance between interpolation points */


    gsl_set_error_handler_off() ;
    Wendland_result result ;
    /* calculating the covariance fct. in the interpolation points 
     * */
    for ( int i=0 ; i < n ; ++i ) {

        points[i] = i*interval ;
        wendland( &result,  points[i]/rnge, mu, smoothness, abstol, reltol ) ;

        if ( check_wendland_errors( &result ) ) {

            wendl[i]    = sill * result.result ;     
        } else {

            return R_NilValue ;
        }
    }

    /* initialisation for the interpolation */
    gsl_interp *interpol = gsl_interp_alloc(gsl_interp_cspline , n ) ;
    gsl_interp_init( interpol,  points, wendl, n ) ;
    gsl_interp_accel *acc =  gsl_interp_accel_alloc() ;

    /* calculating the covariance matrix */
    for ( int i=0 ; i<length ; i++ ) {

        if ( *(p_dist+i) < eps ) {

            REAL(RESULT)[i] = sill + nugget ;
        } else if ( *(p_dist+i) < rnge ) {

            REAL(RESULT)[i] = 
                gsl_interp_eval( interpol, points, wendl, 
                        *(p_dist + i), acc )  ;
        } else {

            REAL(RESULT)[i] = 0 ;
        }
    }

    gsl_interp_accel_free( acc ) ;
    gsl_interp_free( interpol ) ;   
    free(points) ;
    free(wendl) ;
    UNPROTECT(1) ; /* RESULT */

    return RESULT ;
}
