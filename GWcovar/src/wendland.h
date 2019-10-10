/* This file is part of the R-package 'GWcovar'
 *
 * Copyright (C) 2019 Josef Stocker <josef@josefstocker.ch>
 * 
 * 'GWcovar' is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 * */

#ifndef WENDLAND_H_
#define WENDLAND_H_


/* ****************************************************************************
 * ** Include directives  *****************************************************
 * ***************************************************************************/

#include "stddef.h" /* for type size_t */



/* ***************************************************************************
 * ** Public data structures *************************************************
 * **************************************************************************/

typedef struct {
/* ***************************************************************************
 * Struct used to return the value calculated by the function in wendland.c to
 * the program that called the function.
 * **************************************************************************/
    double result ;
    /* result of the numerical integration */

    double abserr ;
    /* error of the numerical integration */
    
    size_t neval ;
    /* number of evaluation for num. integration */

    int error ;
    /* exit status of the GSL function used for numerical integration */
    
    int error_b ;
    /* exit status of the GSL function used for the beta function */
} Wendland_result ;




/* ***************************************************************************
 * ***************************************************************************
 * ** Public functions  ******************************************************
 * ***************************************************************************
 * **************************************************************************/

int 
check_wendland_errors ( 
/* ***************************************************************************
 * The function 'void check_wendland_errors(...)' checks if there has been an
 * errror in the calculation done by a function from 'wendland.c'
 *
 * The function takes an argument 'result' of type 'Wendland_result' and
 * returns '1' if both 'result.error' and 'result.error_b' are '0'. Otherwise
 * it returns '0'. If there is an error, 'check_wendland_errors' also produces
 * an error message using REprintf. 
 * **************************************************************************/
        Wendland_result* result 
        )  ;


void 
wendland(   
/* ***************************************************************************
 * The function 'void wendland(...)' calculates the value of the GW
 * correlation function using non adaptive Gauss-Konrod integration. For the
 * integration the procedure 'int gsl_integration_qng(...)' from the 'GNU
 * Scientific Library' is used.
 *
 *
 *
 * ****************
 * ** Arguments: ** 
 * ****************
 *
 *  ->  Wendland_result* result:    is used to return the resulting value to
 *                          the function from where 'wendland(...)' was called
 * 
 *  ->  double dist:        distance between the two locations for which the GW
 *                          covariance function is calculated
 *
 *  ->  double mu:          parameter of the GW covariance function
 *
 *  ->  double smoothness:  smoothness parameter of the GW covariance function
 *
 *  ->  double abstol:      absolute tolerance for the numerical integration
 *
 *  ->  double reltol:      relative tolerance for the numerical integration
 *
 *
 * ***************************************************************************/
        Wendland_result* result ,
        double dist,
        double mu,
        double smoothness,
        double abstol,
        double reltol
        ) ;


void 
wendland_qag (       
/* ***************************************************************************
 * The function 'double wendland(...)' returns the value of the GW
 * correlation function using adaptive Gauss-Konrod integration. For the
 * integration the procedure 'int gsl_integration_qag(...)' from the 'GNU
 * Scientific Library' is used.
 *
 *
 *
 * **************
 * **Arguments:** 
 * **************
 * 
 * 
 *  ->  Wendland_result* result:    is used to return the resulting value to
 *                          the function from where 'wendland(...)' was called
 *
 *  ->  double dist:        distance between the two locations for 
 *                          which the GW covariance function is 
 *                          calculated
 *
 *  ->  double mu:          parameter of the GW covariance function
 *
 *  ->  double smoothness:  smoothness parameter of the GW covariance
 *                          function
 *
 *  ->  double abstol:      absolute tolerance for the numerical 
 *                          integration
 *
 *  ->  double reltol:      relative tolerance for the numerical 
 *                          integration
 *
 *  ->  int intervals:      the number of intervalls used for the adaptive
 *                          Gauss-Konrod integration
 *
 *  ->  int key:            parameter for the Gauss-Konrod integration
 *
 * ***************************************************************************/
        Wendland_result* result ,
        double dist,
        double mu,
        double smoothness,
        double abstol,
        double reltol,
        int intervals,
        int key
        ) ;

#endif  /* #ifndef WENDLAND_H_ */
