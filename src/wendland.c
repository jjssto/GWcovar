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


/* ***************************************************************************
 * ** Include directives  ****************************************************
 * **************************************************************************/

#include "wendland.h"

#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "stdio.h"
#include "math.h" 
#include "R_ext/Print.h"

/* ***************************************************************************
 * ** Private data structures ************************************************
 * **************************************************************************/
typedef struct {
    /* struct which is used to pass the parameters to the function, that
     * is integrated to get the values of the Wendland covariance
     * function */
    double dist ;
    double mu ;
    double smoothness ;
} Fct_params ;



/* ***************************************************************************
 * ** Functions **************************************************************
 * **************************************************************************/


/* ***********************
 * ** private functions **
 * **********************/

double 
fct1 (    

        double u , 
        void *p 
     ) 
/* function that is integrated to get the GW correlation
 * function */
{
    Fct_params *params = (Fct_params * )p ;
    double dist = (params->dist) ;
    double smoothness = (params->smoothness) ;
    double mu = (params->mu) ;
    return u * pow( (u*u - dist*dist) , smoothness-1.0 ) * pow( 1.0 - u , mu ) ; 
}

double 
fct2 (   

        double u , 
        void *p 
     ) 
/* alternative version of the function, that is integrated
 * to get GW correlation function */
{
    Fct_params *params = (Fct_params * )p ;
    double dist = (params->dist) ;
    double smoothness = (params->smoothness) ;
    double mu = (params->mu) ;
    return pow( (u*u - dist*dist) , smoothness ) * pow( 1.0 - u , mu -1.0 ) ; 
}




/* **********************
 * ** public functions **
 * *********************/
int 
check_wendland_errors ( 
        Wendland_result* result 
        ) 
{
    int ret = 1 ;
    if ( result->error != 0 ) {

        ret = 0 ;
        REprintf(
                "%s\n%s%s\n",
                "Error occured during numerical integration",
                "GSL Error: ",
                gsl_strerror(result->error) 
               ) ;

    }
    if ( result->error_b != 0 ) {

        ret = 0 ;
        REprintf(
                "%s\n%s%s\n",
                "Error occured while calculating the beta function",
                "GSL Error: ",
                gsl_strerror(result->error_b) 
               ) ;

    }
    return ret ;
}

void
wendland (   

        Wendland_result* result ,
        double dist,        /*distance between locations*/
        double mu,          /*param. of correlation function*/
        double smoothness,  /*param. of correlation function*/
        double abstol,      /*param. for integration*/
        double reltol       /*param. for integration*/
        ) 
/* The function 'double wendland([...])' returns the value of the GW
 * correlation function using non adaptive Gauss-Konrod integration. For the
 * integration the procedure 'int gsl_integration_qng([...])' from the 'GNU
 * Scientific Library' is used.
 * */
{
    if ( smoothness == 0) {

        if ( dist < 1 ) {

            result->result = pow( ( (double) 1-dist) , mu )  ;
        } else {

            result->result = 0  ;
        }
        result->error = 0 ;
        result->error_b = 0 ;
    } else {

        if (dist < 1) {

            gsl_set_error_handler_off() ;
            gsl_function F ;
            Fct_params params = { dist, mu, smoothness } ;
            F.function = &fct2 ;
            F.params = &params ;

            result->error = gsl_integration_qng(

                    &F , 
                    dist, // a
                    1.0, // b
                    abstol, //epsabs
                    reltol, //epsrel
                    &(result->result) ,
                    &(result->abserr) , 
                    &(result->neval)
                    ) ;

        } else {

            result->result = 0 ;
            result->error = 0 ;
            result->error_b = 0 ;
        }
        if ( result->result != 0 ) {

            gsl_sf_result result_beta ;
            result->error_b = gsl_sf_beta_e(1+ 2.0*smoothness,  mu, &result_beta ) ;

            if ( result->error_b == 0 ) {

                result->result /= result_beta.val ;
            }
        }
    }
}

void 
wendland_qag (       

        Wendland_result* result ,
        double dist,        /*distance between locations*/
        double mu,          /*param. of correlation func.*/
        double smoothness,  /*param. of correlation func.*/
        double abstol,      /*abs. tolerance for integration*/
        double reltol,      /*rel. tolerance for integration*/
        int intervals,      /*parameter for integration*/
        int key             /*parameter for integration*/
        ) 
/* The function 'double wendland([...])' returns the value of the GW
 * correlation function using adaptive Gauss-Konrod integration. For the
 * integration the procedure 'int gsl_integration_qag([...])' from the 'GNU
 * Scientific Library' is used.
 * */
{
    if ( smoothness == 0) {

        if ( dist < 1 ) {

            result->result =  pow( ( (double) 1-dist) , mu )  ;
        } else {

            result->result = 0  ;
        }
        result->error = 0 ;
        result->error_b = 0 ;
    } else {

        if (dist < 1) {

            gsl_set_error_handler_off();
            gsl_function F ;
            Fct_params params = { dist, mu, smoothness } ;
            F.function = &fct2 ;
            F.params = &params ;

            gsl_integration_workspace *p_workspace = 
                gsl_integration_workspace_alloc( intervals ) ;
            result->error = gsl_integration_qag(

                    &F , 
                    dist, // a
                    1.0, // b
                    abstol, //epsabs
                    reltol, //epsrel
                    intervals,
                    key,
                    p_workspace,
                    &(result->result) ,
                    &(result->abserr) 
                    ) ;
            gsl_integration_workspace_free( p_workspace ) ;


        } else {

            result->result = 0 ;
            result->error = 0 ;
            result->error_b = 0 ;
        }
        if ( result->result != 0 ) {

            gsl_sf_result result_beta ;
            result->error_b = gsl_sf_beta_e(1+ 2.0*smoothness,  mu, &result_beta ) ;

            if ( result->error_b == 0 ) {

                result->result /= result_beta.val ;
            } 
        }    
    }
}
