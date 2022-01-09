
/* 
* separray.h - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/

#ifndef __SEPARRAY_H__
#define __SEPARRAY_H__

/* For error messages and some standard headers */
#include "sepmisc.h"

typedef double*** septensor;
typedef double** sepmatrix;
typedef double* sepvector;


/**
 * Allocates memory to vector array  
 * @param length Length of vector
 * @return Pointer of type DOUBLE to array 
*/
double *sep_vector(size_t length);


/**
 * Allocates memory to matrix array
 * @param nrow Number of rows
 * @param ncol Number of columns  
 * @return Pointer of type DOUBLE to array
 */
double **sep_matrix(size_t nrow, size_t ncol);


/**
 * Allocates memory for three dimensional array.
 * @param n1 Number of elements in dimension one
 * @param n2 Number of elements in dimension two
 * @param n3 Number of elements in dimension three
 * @return Pointer of type DOUBLE to array. 
 */
double ***sep_tensor(size_t n1, size_t n2, size_t n3);


/**
 * Allocates memory to vector array
 * @param length Length of vector
 * @return Pointer of type INT to array
 */
int *sep_vector_int(size_t length);

/**
 * Allocates memory to vector array
 * @param length Length of vector
 * @return Pointer of type SIZE_T to array
 */
size_t *sep_vector_size_t(size_t length);

/**
 * Allocates memory to matrix array
 * @param nrow Number of rows
 * @param ncol Number of columns
 * @return Pointer of type INT to array 
 */
int **sep_matrix_int(size_t nrow, size_t ncol);


/**
 * Allocates memory for three dimensional array.
 * @param n1 Number of elements in dimension one
 * @param n2 Number of elements in dimension two
 * @param n3 Number of elements in dimension three
 * @return Pointer of type INT to array. 
 */
int ***sep_tensor_int(size_t n1, size_t n2, size_t n3);


/**
 * Frees the memory of type DOUBLE
 * @param **ptr The pointer pointing to memory
 * @param nrow Number of rows 
 */
void sep_free_matrix(double **ptr, size_t nrow);


/**
 * Frees the memory allocated be sep_tensor().
 * @param ***a Pointer to array
 * @param n1  Number of elements in dimension one
 * @param n2  Number of elements in dimension two
 */
void sep_free_tensor(double ***a, size_t n1, size_t n2);


/**
 * Frees the memory of type INT
 * @param **ptr The pointer pointing to memory
 * @param nrow Number of rows 
 */
void sep_free_matrix_int(int **ptr, size_t nrow);


/**
 * Frees the memory of type INT
 * @param ***ptr The pointer pointing to memory
 * @param n1 Number of elements dimension one
 * @param n2 Number of elements dimension two
 */
void sep_free_tensor_int(int ***ptr, size_t n1, size_t n_2);

/**
 * Writes vector to file
 * @param *vec Pointer to vector of type DOUBLE
 * @param file[] Name of file
 * @param opt Option: "w" for writing "a" for appending
 * @param length Length of vector
 */
void sep_fprintf_vector(double *vec, const char file[], char opt[2], size_t length);

/**
 * Writes matrix to file
 * @param **A Pointer to matrix array
 * @param file[] Name of file
 * @param opt Option: "w" for writing "a" for appending
 * @param nrow Number of rows of matrix
 * @param ncol Number of columns of matrix
 */
void sep_fprintf_matrix(double **A, const char file[], char opt[2], 
			size_t nrow, size_t ncol);


/**
 * Writes vector to file
 * @param *vec Pointer to vector of type DOUBLE
 * @param file[] Name of file
 * @param opt Option: "w" for writing "a" for appending
 * @param length Length of vector
 */
void sep_fprintf_vector_int(int *vec, const char file[], char opt[2], size_t length);

/**
 * Writes matrix to file
 * @param **A Pointer to matrix array
 * @param file[] Name of file
 * @param opt Option: "w" for writing "a" for appending
 * @param nrow Number of rows of matrix
 * @param ncol Number of columns of matrix
 */
void sep_fprintf_matrix_int(int **A, const char file[], char opt[2], 
			    size_t nrow, size_t ncol);

/**
 * Writes vector to stdout
 * @param *a Pointer to vector of type DOUBLE
 * @param length Length of vector
 */
void sep_printf_vector(double *a, int length);


/**
 * Writes matrix to stdout
 * @param **a Pointer to matrix array
 * @param nrow Number of rows of matrix
 * @param ncol Number of columns of matrix
 */
void sep_printf_matrix(double **a, size_t nrow, size_t ncol);

/**
 * Writes vector to stdout
 * @param *a Pointer to vector of type int
 * @param length Length of vector
 */
void sep_printf_vector_int(int *a, int length);

/**
 * Writes matrix to stdout
 * @param **a Pointer to matrix array
 * @param nrow Number of rows of matrix
 * @param ncol Number of columns of matrix
 */
void sep_printf_matrix_int(int **a, size_t nrow, size_t ncol);


/**
 * Sets all elements in vector to certain value
 * @param *a Pointer to vector array
 * @param length Length of array
 * @param Value of elements
 */
void sep_vector_set(double *a, size_t length, double value);

/** 
 * Sets all elements in matrix to certain value
 * @param **a Pointer to matrix array
 * @param nrow Number of rows in matrix
 * @param ncol Number of columns in matrix
 * @param Value of elements
 */
void sep_matrix_set(double **a, size_t nrow, size_t ncol, double value);

/** Sets all elements in tensor to certain value
 * @param ***a Pointer to tensor array
 * @param n1 Number of elements in dimension one
 * @param n2 Number of elements in dimension two
 * @param n3 Number of elements in dimension three
 * @param Value of elements
 */
void sep_tensor_set(double ***a, size_t n1, size_t n2, size_t n3,
		    double value);

/**
 * Sets all elements in vector to certain value
 * @param *a Pointer to vector array
 * @param length Length of array
 * @param Value of elements
 */
void sep_vector_int_set(int *a, size_t length, int value);

/** 
 * Sets all elements in matrix to certain value
 * @param **a Pointer to matrix array
 * @param nrow Number of rows in matrix
 * @param ncol Number of columns in matrix
 * @param Value of elements
 */
void sep_matrix_int_set(int **a, size_t nrow, size_t ncol, int value);


/** Sets all elements in tensor to certain value
 * @param ***a Pointer to tensor array
 * @param n1 Number of elements in dimension one
 * @param n2 Number of elements in dimension two
 * @param n3 Number of elements in dimension three
 * @param Value of elements
 */
void sep_tensor_int_set(int ***a, size_t n1, size_t n2, size_t n3,
			int value);

#ifdef COMPLEX			
typedef _Complex double*** sepctensor;
typedef _Complex double** sepcmatrix;
typedef _Complex double* sepcvector;


/**
 * Allocates memory to COMPLEX DOUBLE vector array  
 * @param length Length of vector
 * @return Pointer of type COMPLEX DOUBLE to array 
*/
_Complex double *sep_complex_vector(size_t length);

/**
 * Allocates memory to COMPLEX DOUBLE matrix array
 * @param nrow Number of rows
 * @param ncol Number of columns  
 * @return Pointer of type COMPLEX DOUBLE to array
 */
_Complex double **sep_complex_matrix(size_t nr, size_t nc);


/**
 * Allocates memory for three dimensional COMPLEX DOUBLE array.
 * @param n1 Number of elements in dimension one
 * @param n2 Number of elements in dimension two
 * @param n3 Number of elements in dimension three
 * @return Pointer of type COMPLEX DOUBLE to array. 
 */
_Complex double ***sep_complex_tensor(size_t n1, size_t n2, size_t n3);


/**
 * Frees the memory of type COMPLEX DOUBLE
 * @param a triple pointer to array
 * @param n1 First dimension length
 * @param n2 Second dimension length
 */
void sep_free_complex_tensor(_Complex double ***a, size_t n1, size_t n2);

/**
 * Frees the memory of type COMPLEX DOUBLE
 * @param a The pointer pointing to memory
 * @param nr Number of rows 
 */
void sep_free_complex_matrix(_Complex double **a, size_t nr);

#ifndef DOXYGEN_SKIP
void sep_fprintf_complex_vector(_Complex double *a, const char *file, 
				char opt[2], size_t length);
void sep_fprintf_complex_matrix(_Complex double **a, const char *file, 
				char opt[2], size_t nr, size_t nc);
#endif
#endif

#ifndef DOXYGEN_SKIP

void sep_free_tensor_float(float ***a, size_t n1, size_t n2);
void sep_free_matrix_float(float **ptr, size_t nrow);
float ***sep_tensor_float(size_t n1, size_t n2, size_t n3);
float **sep_matrix_float(size_t nrow, size_t ncol);
								       
#endif
#endif

