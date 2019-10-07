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


/* ****************************************************************************
 * ** Errors ******************************************************************
 * ***************************************************************************/
#define INVALID_ARG 50
#define NULL_PTR 51



/* ***************************************************************************
 * ** Public data structures *************************************************
 * **************************************************************************/
typedef struct {
/* ***************************************************************************
 * ...
 * **************************************************************************/
    double result ;
    /* result of integration */

    double abserr ;
    /* result of integration */
    
    size_t neval ;
    /* number of integration */

    int error ;
    /* 0 if no error occured, otherwise GSL error code */
    
    int error_b ;
    /* 0 if no error occured, otherwise GSL error code */
} Wendland_result ;




/* ***************************************************************************
 * ***************************************************************************
 * ** Public functions  ******************************************************
 * ***************************************************************************
 * **************************************************************************/

int 
check_wendland_errors ( 
/* ***************************************************************************
 * ...
 * **************************************************************************/
        Wendland_result* result 
        )  ;


void 
wendland(   
/* ***************************************************************************
 * The function 'double wendland([...])' returns the value of the GW
 * correlation function using non adaptive Gauss-Konrod integration. For the
 * integration the procedure 'int gsl_integration_qng([...])' from the 'GNU
 * Scientific Library' is used.
 *
 *
 *
 * ****************
 * ** Arguments: ** 
 * ****************
 * 
 *  ->  double dist:        distance between the two locations for 
 *                          which the GW covariance function is 
 *                          calculated
 *
 *  ->  double mu:          parameter of the GW covariance 
 *                          function
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
 * The function 'double wendland([...])' returns the value of the GW
 * correlation function using adaptive Gauss-Konrod integration. For the
 * integration the procedure 'int gsl_integration_qag([...])' from the 'GNU
 * Scientific Library' is used.
 *
 *
 *
 * **************
 * **Arguments:** 
 * **************
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
