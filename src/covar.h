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

#ifndef COVAR_H_
#define COVAR_H_


/* ***************************************************************************
 * ** Errors *****************************************************************
 * **************************************************************************/

#define INVALID_PARAMETERS 60
#define METHOD_NOT_DEFINED 61
#define WRONG_QAG_PARAMETERS 62
#define WRONG_INTERPOL_PARAMETER 63
#define INVALID_TOLERANCE 64


/* ***************************************************************************
 * ** Includes  **************************************************************
 * **************************************************************************/

#include "Rinternals.h"



/* ***************************************************************************
 * ** Functions  *************************************************************
 * **************************************************************************/
SEXP 
covar_m_dist (  
/* **************************************************************************** 
 * The function 'SEXP covar_m_dist(...)' returns the Generalized Wendland
 * (GW) covariance matrix based on a distance matrix in standard R matrix
 * format.  The values of the covariance function are calculated using  
 * non-addaptive Gauss-Kronrod integration.
 *
 * 
 *  ****************
 *  ** Arguments: **
 *  ****************
 *  
 *  ->  SEXP DIST:      Distance matrix in standard R matrix format. If 'DIST'
 *                      is a square matrix, only the upper triangular matrix 
 *                      is accessed.  
 *
 *  -> SEXP MU:         Parameter of the GW covariance function
 *
 *  -> SEXP SMOOTHNESS: Parameter of the GW covariance function
 *
 *  -> SEXP SILL:       Parameter of the GW covariance function. 'SILL'
 *                      controls the variance at the locations
 *
 *  -> SEXP RNGE:       Parameter of the GW covariance function. 'RNGE'
 *                      controls the radius of the compact support of the GW
 *                      covariance funcion.
 *
 *  -> SEXP NUGGET:     Parameter of the GW covariance function. 'NUGGET'
 *                      controls the nugget of the GW covariance function.
 *
 *  -> SEXP ABSTOL:     Parameter for the numerical integration: absolute
 *                      tolerance.
 *
 *  -> SEXP RELTOL:     Parameter for the numerical integration: relative
 *                      tolerance.
 *
 *  ->  SEXP NBR_INTERPOL:  Number of equidistant points in which the GW
 *                      covariance function is evaluated, if interpolation is
 *                      used.
 *                          
 * ****************************************************************************/
        SEXP DIST ,      /* distance matrix */
        SEXP MU ,        /* param. of GW correlation function */
        SEXP SMOOTHNESS ,/* param. of GW correlation function */
        SEXP SILL ,      /* param. of GW correlation function */
        SEXP RNGE ,      /* param. of GW correlation function */
        SEXP NUGGET,     /* param. of GW correlation function */
        SEXP ABSTOL,     /* absolute tolerance for integration */
        SEXP RELTOL      /* relative tolerance for integration */
        ) ;

SEXP 
covar_interpol (
/* **************************************************************************** 
 * The function 'SEXP covar_interpol(...)' returns the Generalized Wendland
 * (GW) covariance matrix based on a distance matrix in standard R matrix
 * format.  The values of the covariance function are interpolated using cubic
 * splines. For the calculation of the GW covariance function the
 * non-addaptive Gauss-Kronrod algorithm is used.
 *
 * 
 *  ****************
 *  ** Arguments: **
 *  ****************
 *  
 *  ->  SEXP DIST:      Distance matrix in standard R matrix format. If 'DIST'
 *                      is a square matrix, only the upper triangular matrix 
 *                      is accessed.  
 *
 *  -> SEXP MU:         Parameter of the GW covariance function
 *
 *  -> SEXP SMOOTHNESS: Parameter of the GW covariance function
 *
 *  -> SEXP SILL:       Parameter of the GW covariance function. 'SILL'
 *                      controls the variance at the locations
 *
 *  -> SEXP RNGE:       Parameter of the GW covariance function. 'RNGE'
 *                      controls the radius of the compact support of the GW
 *                      covariance funcion.
 *
 *  -> SEXP NUGGET:     Parameter of the GW covariance function. 'NUGGET'
 *                      controls the nugget of the GW covariance function.
 *
 *  -> SEXP ABSTOL:     Parameter for the numerical integration: absolute
 *                      tolerance.
 *
 *  -> SEXP RELTOL:     Parameter for the numerical integration: relative
 *                      tolerance.
 *
 *  ->  SEXP NBR_INTERPOL:  Number of equidistant points in which the GW
 *                      covariance function is evaluated, if interpolation is
 *                      used.
 *
 * ****************************************************************************/
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
        ) ;


SEXP covar_vector_dir (
/* ****************************************************************************
 * The function 'int covar_vector_dir  (...)' calculates the Generalized
 * Wendland (GW) covariance matrix by modifying the R vector 'DIST'. The
 * function returns 0 if the calculation finished without erros.  Depending on
 * the parameters used the integral is calculated with the non-adaptive
 * Gauss-Kronrod algorithm using the 'GNU Scientific Library'. 
 *
 * 
 *  ****************
 *  ** Arguments: **
 *  ****************
 *  
 *  -> SEXP DIST:       Distance matrix in standard R matrix format. If 'DIST'
 *                      is a square matrix, only the upper triangular matrix 
 *                      is accessed.  
 *
 *  -> SEXP LENGTH:     The number of elements in DIST.
 *
 *  -> SEXP MU:         Parameter of the GW covariance function
 *
 *  -> SEXP SMOOTHNESS: Parameter of the GW covariance function
 *
 *  -> SEXP SILL:       Parameter of the GW covariance function. 'SILL'
 *                      controls the variance at the locations
 *
 *  -> SEXP RNGE:       Parameter of the GW covariance function. 'RNGE'
 *                      controls the radius of the compact support of the GW
 *                      covariance funcion.
 *
 *  -> SEXP NUGGET:     Parameter of the GW covariance function. 'NUGGET'
 *                      controls the nugget of the GW covariance function.
 *
 *  -> SEXP ABSTOL:     Parameter for the numerical integration: absolute
 *                      tolerance.
 *
 *  -> SEXP RELTOL:     Parameter for the numerical integration: relative
 *                      tolerance.
 *
 *  -> SEXP EPS:        Treshold below which a number is considered to be equal
 *                      to zero.
 *
 * ****************************************************************************/
        SEXP DIST ,         /* R vector containing distances */    
        SEXP LENGTH ,       /* length of 'SEXP DIST' */ 
        SEXP MU ,           /* param. of the GW covariance fct */
        SEXP SMOOTHNESS ,   /* param. of the GW covariance fct */
        SEXP SILL ,         /* param. of the GW covariance fct */
        SEXP RNGE ,         /* param. of the GW covariance fct */
        SEXP NUGGET ,       /* param. of the GW covariance fct */
        SEXP ABSTOL ,       /* abs. tolerance for integration */
        SEXP RELTOL ,       /* rel. tolerance for integration */
        SEXP EPS           /* treshhold below which values are */
        ) ;

SEXP covar_vector_interpol (
/* *****************************************************************************
 * The function 'int covar_vector_dir  (...)' calculates the Generalized Wendland
 * (GW) covariance matrix by modifying the R vector 'DIST'. The function
 * returns 0 if the calculation finished without erros.  Depending on the
 * parameters used the integral is calculated with the non-adaptive
 * Gauss-Kronrod algorithm using the 'GNU Scientific Library'. 
 *
 * 
 *  ****************
 *  ** Arguments: **
 *  ****************
 *  
 *  -> SEXP DIST:       Distance matrix in standard R matrix format. If 'DIST'
 *                      is a square matrix, only the upper triangular matrix 
 *                      is accessed.  
 *
 *  -> SEXP LENGTH:     The number of elements in DIST.
 *
 *  -> SEXP MU:         Parameter of the GW covariance function
 *
 *  -> SEXP SMOOTHNESS: Parameter of the GW covariance function
 *
 *  -> SEXP SILL:       Parameter of the GW covariance function. 'SILL'
 *                      controls the variance at the locations
 *
 *  -> SEXP RNGE:       Parameter of the GW covariance function. 'RNGE'
 *                      controls the radius of the compact support of the GW
 *                      covariance funcion.
 *
 *  -> SEXP NUGGET:     Parameter of the GW covariance function. 'NUGGET'
 *                      controls the nugget of the GW covariance function.
 *
 *  -> SEXP ABSTOL:     Parameter for the numerical integration: absolute
 *                      tolerance.
 *
 *  -> SEXP RELTOL:     Parameter for the numerical integration: relative
 *                      tolerance.
 *
 *  -> SEXP EPS:        Treshold below which a number is considered to be equal
 *                      to zero.
 *
 * 
 *  ->  SEXP NBR_INTERPOL:  Number of equidistant points in which the GW
 *                      covariance function is evaluated, if interpolation is
 *                      used.
 *
 * ****************************************************************************/
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
        ) ;

#endif  /* COVAR_H_ */
