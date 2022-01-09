/* 
* seputil.h - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/

#ifndef __SEPUTIL_H__
#define __SEPUTIL_H__

#include "sepmisc.h"
#include "separray.h"

/**
 *  Solves a set of linear equations of the form Ax = b
 *  using LU decomposition  - taken from Rapaport
 *  @param b is the known n length vector
 *  @param A is a n times n matrix 
 *  @param n is the length number of unkowns.
 *  Notice: b and A will be overwritten. The function returns 
 *  the result in b
 *  
 */
void sep_solvelineq(double *b, double **A, int n);


/**
 *  Solves a set of linear equations of the form Ax = b
 *  using LU decomposition  - taken from Rapaport
 *  @param ret the return vector
 *  @param b is the known n length vector
 *  @param A is a n times n matrix 
 *  @param n is number of unkowns.
 *  
 */
void sep_solvelineq1(double *ret, double *x, double **A, int n);

/**
 * Fits linear model to data: data = a0*f0(x) + a2*f2(x)+....aN*fN(x)
 * @param param Array that holds the returning fitted parameters
 * @param x  Array with the independent variables 
 * @param data Data array to be fitted 
 * @param ndata Number of data point (same as length of x) 
 * @param nparam Number of parameters minus 1 (N-1 in the above equation).  
 * @param model Function descriping the model to be fitted without coeffcients. 
 * The first argument is the value of the
 * independent variable (x) the second is the term number. E.g. the
 * model is: f(x) = a0*x^2 + a1*cos(x) the the function is 
 * <br> 
 * double fun(double x, int n){ <br>
 * y = 0.0; <br>
 * if (n==0) y = x*x; <br>
 * else if ( n==1 )
 * y = cos(x); <br>
 * return y; 
 * <br>}<br>
 */
void sep_linfit(double *param, double *x, double *data,  
		size_t ndata,  size_t nparam, 
		double (*model)(double, int));

/**
 * Same as sep_linfit but where the terms are given in a lists.
 * @param param Is a vector that holds the returning fitted 
 * parameters 
 * @param x Is a vector holding the independent variables
 * @param data Data array to be fitted
 * @param ndata Number of data point (same as length of x)
 * @param *termlist Array holding the terms e.g. [0,2,10]
 * @param nterms Length of termlist array.
 * @param model Function descriping the model to be fitted without 
 * coeffcients. Also, see sep_linfit() 
.*/
void sep_linfit_list(double *param, double *x, 
		     double *data, size_t ndata,  
		     unsigned int *termlist, size_t nterm,
		     double (*fun)(double, int));

/** 
 * Calculates the mean of a data array
 * @param *vec Pointer to array
 * @param length Length of array
 * @return Mean of array
 */
double sep_vector_mean(double *vec, int length);

/** 
 * Calculates the standard deviation of a data array
 * @param *vec Pointer to array
 * @param length Length of array
 * @return Standard deviation of array
 */
double sep_vector_std(double *vec, int length);

/** 
 * Calculates the skewness of a data array
 * @param *vec Pointer to array
 * @param length Length of array
 * @return Skewness of array
 */
double sep_vector_skewness(double *vec, int length);

/** 
 * Calculates the kurtosis of a data array
 * @param *vec Pointer to array
 * @param length Length of array
 * @return Kurtosis of array
 */
double sep_vector_kurtosis(double *vec, int length);


/**
 * Calculates the linear correlation coefficient between two data
 * sets x and y. 
 * @param *x First data array
 * @param *y Second data array
 * @param ndata Length of data arrays
 * @return The linear correlation coefficient
 */
double sep_lincorr(double *x, double *y, size_t ndata);

/** 
 * Finding the maximum of a data array
 * @param *vec Pointer to array
 * @param length Length of array
 * @return the maximum of array
 */
double sep_vector_max(double *vec, int length);

/** 
 * Finding the minimum of a data array
 * @param *vec Pointer to array
 * @param length Length of array
 * @return the minimum of array
 */
double sep_vector_min(double *vec, int length);

/**
 * Calculates the autoscorrelation of a data array
 * @param *out Pointer to autocorrelation data array
 * @param *in Pointer to data array
 * @param lvec Length of array
 */
void sep_autocorr(double *out, double *in, size_t lvec);

/** 
 * Calculates the integral of the data array from a to b using
 * the trapziodal rule
 * @param *f Pointer to array
 * @param a Lower limit
 * @param b Upper limit
 * @param length Length of array
 */
double sep_trapz(double *f, double a, double b, int length);


/**
 * The sign function 
 * @param x Argument
 * @return the sign of argument
 */
 int sep_sgn(double x);


/** 
 * The Heaviside function 
 * @param x Argument
 * @return Zero if x is negative else one.
 */
 int sep_heaviside(double x);

/**
 * Calculates the dot product of two vectors
 * @param *a First vector 
 * @param *b Second vector
 * @param length Length of vectors
 * @return Dot product. 
*/
double sep_dot(double *a, double *b, int length);

/**
* Calculates the norm of a vector
* @param *a Vector
* @param length Length of vector
* @return The norm 
*/
 double sep_vector_norm(double *a, int length);

/**
* Calculates the trace of a square matrix
* @param **A Matrix
* @param nrow Number of rows (and columns)
* @return The trace
*/  
 double sep_trace(double **A, size_t nrow);

/**
 * Calculates the vector sum
 * @param a The vector
 * @param length Number of vector elements
 * @return  The sum
 */
 double sep_vector_sum(double *a, size_t length);

/**
 * Multiply a vector with a scalar elementwise. The vector values are 
 * overwritten
 * @param a The vector
 * @param value The scalar
 * @param lvec Number of vector elements
 */
 void sep_vector_multiply(double *a, double value, size_t lvec);


/**
 * Add a vector with a scalar elementwise. The vector values are 
 * overwritten
 * @param a The vector
 * @param value The scalar
 * @param lvec Number of vector elements
 */
 void sep_vector_add(double *a, double value, size_t lvec);


/**
 * Multiply a matrix with a scalar elementwise. The matrix values are 
 * overwritten
 * @param A The matrix pointer
 * @param value The scalar
 * @param nrow Number of rows
 * @param ncol Number of columns
 */
 void sep_matrix_multiply(double **A, double value, 
                                size_t nrow, size_t ncol);


/**
 * Add a matrix with a scalar elementwise. The matrix values are 
 * overwritten
 * @param A The matrix pointer
 * @param value The scalar
 * @param nrow Number of rows
 * @param ncol Number of columns
 */
 void sep_matrix_add(double **A, double value, 
                           size_t nrow, size_t ncol);

/** 
 * Doing the standard 3-D cross product
 * @param ret Return array (length 3)
 * @param a 'Left' vector array (length 3)
 * @param b 'Right' vector array (length 3) 
 */
void sep_cross3(double *ret, double *a, double *b);

/**
 * Standard determinant for a real 3 by 3 matrix
 * @param A Real 3 by 3 matrix
 * @return The determinant
 */
double sep_det3(double **A);

/**
 * Fill vector array with values from [a;b] evenly and linearly spaced
 * @param vec Return vector array of length lvec
 * @param a Lower limit
 * @param b Upper limit
 * @param lvec Length of vector array
 */
void sep_linspace(double *vec, double a, double b, size_t lvec);

/**
 * Decompose a rank-2 tensor into av. trace, symmetric and anti-symmetric
 * parts.
 * @param p The average trace (return)
 * @param Pos Symmetric tensor part (return)
 * @param Pa Anti-symmetric tensor part (return)
 * @param P The input rank 2 tensor
 * @param nrc Number of rows (same as number of columns)
 */
void sep_decompose_tensor(double *p, double **Pos, double **Pa, double **P, size_t nrc);

/**
 * Straight forward matrix vector multiplication 
 * @param out Resulting vector array
 * @param A Matrix
 * @param x Vector 
 * @param nr Number of matrix rows
 * @param nc Number of matrix columns 
 */
void sep_mv_multi(double *out, double **A, double *x, int nr, int nc); 

/**
 * Calculates the eigenvalues of a real, 3 by 3, symmetric matrix
 * @param eig Return vector array (length 3)
 * @param A Real, symmetric, 3 by 3 matrix
 */
void sep_eig_real_symmetric(double *eig, double** A);

/** 
 * Evaluates the correlation between two complex arrays with full 
 * delay/lag range. Parallelized with omp.
 * @param A Output matrix array
 * @param B Input first matrix array
 * @param C Input second matrix array
 * @param nrows Number of rows
 * @param ncol Number of columns
 */ 
void sep_corr_complex_arrays(_Complex double **A, _Complex double **B,
                             _Complex double **C, size_t nrows, size_t ncol);



#ifndef DOXYGEN_SKIP
/**
 * Power method for estimating the maximum/minimum eigenvalue for a matrix A
 * Taken from Kreyszig, 9th edition, sect 20.8. 
 */ 
unsigned sep_maxeigvalue_pow(double *q, double **A, int nrc, double eps0, 
			     const unsigned countermax);


int sep_compare_int_ascend(const void *a, const void *b);
int sep_binomial(double p, int n);

#endif

#endif
