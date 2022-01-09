
/* 
* separray.c - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/

#include "separray.h"

double *sep_vector(size_t length){
  size_t n;
  double *ptr;

  ptr = malloc(length*sizeof(double));
  if (ptr == NULL)
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,	__LINE__);
  
  for (n=0; n<length; n++)
    ptr[n] = 0.0;

  return ptr;
}


double **sep_matrix(size_t nrow, size_t ncol){
  double **ptr;
  size_t n, m;

  ptr = malloc(nrow*sizeof(double *));
  if (ptr == NULL)
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
	      __LINE__);

  for (n=0; n<nrow; n++){
    ptr[n] = malloc(ncol*sizeof(double));
    if (ptr[n] == NULL)
      sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
		__LINE__);
  }

  for (n=0; n<nrow; n++)
    for (m=0; m<ncol; m++)
      ptr[n][m] = 0.0;

  return ptr;
}



double ***sep_tensor(size_t n1, size_t n2, size_t n3){
  double ***a;
  size_t n, m, k;

  a = malloc(n1*sizeof(double **));
  if ( a == NULL )
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
	      __LINE__);
    
  for ( n=0; n<n1; n++ ){
    a[n] = sep_matrix(n2, n3);
    if ( a[n] == NULL )
      sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
		__LINE__);
  }

  for ( n=0; n<n1; n++ )
    for ( m=0; m<n2; m++ )
      for ( k=0; k<n3; k++ )
	a[n][m][k] = 0.0;

  return a;
      
}


int *sep_vector_int(size_t length){
  size_t n;
  int *ptr;

  ptr = malloc(length*sizeof(int));
  
  if (ptr == NULL)
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
	      __LINE__);

  for (n=0; n<length; n++)
    ptr[n] = 0;

  return ptr;
}

size_t *sep_vector_size_t(size_t length){
  size_t n;
  size_t *ptr;

  ptr = malloc(length*sizeof(size_t));
  if (ptr == NULL)
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
	      __LINE__);

  for (n=0; n<length; n++)
    ptr[n] = 0;

  return ptr;
}

int **sep_matrix_int(size_t nrow, size_t ncol){
  int **ptr;
  size_t n, m;

  ptr = malloc(nrow*sizeof(int *));
  
  if (ptr == NULL)
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
	      __LINE__);

  for (n=0; n<nrow; n++){
    ptr[n] = malloc(ncol*sizeof(int));
    if (ptr[n] == NULL)
      sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
		__LINE__);
  }

  for (n=0; n<nrow; n++)
    for (m=0; m<ncol; m++)
      ptr[n][m] = 0;

  return ptr;
}



int ***sep_tensor_int(size_t n1, size_t n2, size_t n3){
  int ***a;
  size_t n, m, k;

  a = malloc(n1*sizeof(int **));
  if ( a == NULL )
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
	      __LINE__);
  

  for ( n=0; n<n1; n++ ){
    a[n] = sep_matrix_int(n2, n3);
    if ( a[n] == NULL )
      sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
		__LINE__);
    
  }

  for ( n=0; n<n1; n++ )
    for ( m=0; m<n2; m++ )
      for ( k=0; k<n3; k++ )
	a[n][m][k] = 0;

  return a;
      
}


float **sep_matrix_float(size_t nrow, size_t ncol){
  float **ptr;
  size_t n, m;

  ptr = malloc(nrow*sizeof(float *));
  
  if (ptr == NULL)
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
	      __LINE__);

  for (n=0; n<nrow; n++){
    ptr[n] = malloc(ncol*sizeof(float));
    if (ptr[n] == NULL)
      sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
		__LINE__);
  }

  for (n=0; n<nrow; n++)
    for (m=0; m<ncol; m++)
      ptr[n][m] = 0.0;

  return ptr;
}


float ***sep_tensor_float(size_t n1, size_t n2, size_t n3){
  float ***a;
  size_t n, m, k;

  a = malloc(n1*sizeof(float **));
  if ( a == NULL )
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,__LINE__);

  for ( n=0; n<n1; n++ ){
    a[n] = sep_matrix_float(n2, n3);
    if ( a[n] == NULL )
      sep_error("%s at line %d: Couldn't allocate memory\n", __func__,	__LINE__);
    
  }

  for ( n=0; n<n1; n++ )
    for ( m=0; m<n2; m++ )
      for ( k=0; k<n3; k++ ) a[n][m][k] = 0.0;

  return a;
      
}

void sep_free_matrix(double **ptr, size_t nrow){
  size_t n;

  for (n=0; n<nrow; n++)
    free(ptr[n]);

  free(ptr);
}

void sep_free_tensor(double ***a, size_t n1, size_t n2){
  size_t n;
  
  for ( n=0; n<n1; n++ )
    sep_free_matrix(a[n], n2);

  free(a);  

}


void sep_free_matrix_int(int **ptr, size_t nrow){
  size_t n;

  for (n=0; n<nrow; n++)
    free(ptr[n]);

  free(ptr);
}


void sep_free_tensor_int(int ***a, size_t n1, size_t n2){
  size_t n;
  
  for ( n=0; n<n1; n++ )
    sep_free_matrix_int(a[n], n2);

  free(a);  

}



void sep_free_matrix_float(float **ptr, size_t nrow){
  size_t n;

  for (n=0; n<nrow; n++)
    free(ptr[n]);

  free(ptr);
}


void sep_free_tensor_float(float ***a, size_t n1, size_t n2){
  size_t n;
  
  for ( n=0; n<n1; n++ )
    sep_free_matrix_float(a[n], n2);

  free(a);  

}


void sep_fprintf_vector(double *vec, const char file[], char opt[2], size_t length){
  FILE *fout;
  size_t n;

  fout = fopen(file, opt);
  if (fout == NULL)
    sep_error("%s at line %d: Couldn't open file\n",
	      __func__, __LINE__);
  
  for (n=0; n<length; n++)
    fprintf(fout, "%f ", vec[n]);
  
  fprintf(fout, "\n");

  fclose(fout);
}


void sep_fprintf_vector_int(int *vec, const char file[], 
			    char opt[2], size_t length){
  FILE *fout;
  size_t n;

  fout = fopen(file, opt);
  if (fout == NULL)
    sep_error("%s at line %d: Couldn't open file\n",
	      __func__, __LINE__);
  
  for (n=0; n<length; n++)
    fprintf(fout, "%d ", vec[n]);
  
  fprintf(fout, "\n");

  fclose(fout);
}


void sep_fprintf_matrix(double **A, const char file[], char opt[2], 
			size_t nrow, size_t ncol){
  FILE *fout;
  size_t n, m;

  fout = fopen(file, opt);
  if (fout == NULL)
    sep_error("%s at line %d: Couldn't open file\n",
	      __func__, __LINE__);
      
  for (n=0; n<nrow; n++){
    for (m=0; m<ncol; m++){
      fprintf(fout, "%f ", A[n][m]);
    }
    fprintf(fout,"\n");
  }
  
  fclose(fout);
}


void sep_fprintf_matrix_int(int **A, const char file[], char opt[2], 
			    size_t nrow, size_t ncol){
  FILE *fout;
  size_t n, m;

  fout = fopen(file, opt);
  if (fout == NULL)
    sep_error("%s at line %d: Couldn't open file\n",
	      __func__, __LINE__);
      
  for (n=0; n<nrow; n++){
    for (m=0; m<ncol; m++){
      fprintf(fout, "%d ", A[n][m]);
    }
    fprintf(fout,"\n");
  }
  
  fclose(fout);
}



void sep_printf_matrix(double **a, size_t nrow, size_t ncol){
  size_t n, m;
  
  for ( n=0; n<nrow; n++ ){
    for ( m=0; m<ncol; m++ ){
      printf("%f ", a[n][m]);
    }
    printf("\n");
  }
  
  SEP_FLUSH;
  
}

void sep_printf_vector(double *a, int length){
  int n;

  for ( n=0; n<length; n++ ){
    printf("%f ", a[n]);
  }
  printf("\n");
  
  SEP_FLUSH;
}

void sep_printf_matrix_int(int **a, size_t nrow, size_t ncol){
  size_t n, m;

  for ( n=0; n<nrow; n++ ){
    for ( m=0; m<ncol; m++ ){
      printf("%d ", a[n][m]);
    }
    printf("\n");
  }
  
  SEP_FLUSH;
  
}

void sep_printf_vector_int(int *a, int length){
  int n;

  for ( n=0; n<length; n++ ){
    printf("%d ", a[n]);
  }
  printf("\n");
  
  SEP_FLUSH;
}



void sep_vector_set(double *a, size_t length, double value){
  size_t n;

  for ( n=0; n<length; n++ )
    a[n] = value;

}


void sep_matrix_set(double **a, size_t nrow, size_t ncol, double value){
  size_t n, m;

  for ( n=0; n<nrow; n++ )
    for ( m=0; m<ncol; m++ )
      a[n][m] = value;
}


void sep_tensor_set(double ***a, size_t n1, size_t n2, size_t n3,
		    double value){
  size_t n, m, k;

  for (n =0; n<n1; n++ )
    for ( m=0; m<n2; m++ )
      for ( k=0; k<n3; k++ )
	a[n][m][k] = value;
  
  
}


void sep_vector_int_set(int *a, size_t length, int value){
  size_t n;

  for ( n=0; n<length; n++ )
    a[n] = value;

}


void sep_matrix_int_set(int **a, size_t nrow, size_t ncol, int value){
  size_t n, m;

  for ( n=0; n<nrow; n++ )
    for ( m=0; m<ncol; m++ )
      a[n][m] = value;
}


void sep_tensor_int_set(int ***a, size_t n1, size_t n2, size_t n3,
			int value){
  size_t n, m, k;

  for (n =0; n<n1; n++ )
    for ( m=0; m<n2; m++ )
      for ( k=0; k<n3; k++ )
				
	a[n][m][k] = value;
  
}

#ifdef COMPLEX

_Complex double *sep_complex_vector(size_t length){
		
  _Complex double *retval = malloc(length*sizeof(_Complex double));
  if ( retval == NULL ){
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  }

  for ( size_t i=0; i<length; i++ )
    retval[i] = 0.0 + SEP_I*0.0;
	
  return retval;
}

_Complex double **sep_complex_matrix(size_t nr, size_t nc){
  size_t i, j;
	
  _Complex double **retval = malloc(nr*sizeof(_Complex double)); 
  if ( retval == NULL )
    sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
	
  for ( i=0; i<nr; i++ ){
    retval[i] = malloc(nc*sizeof(_Complex double));
    if ( retval[i] == NULL )
      sep_error("%s at line %d: Couldn't allocate memory", __func__, __LINE__);
  }

  for ( i=0; i<nr; i++ )
    for ( j=0; j<nc; j++ )
      retval[i][j] = 0.0 + 0.0*SEP_I;
	
  return retval;
}


_Complex double ***sep_complex_tensor(size_t n1, size_t n2, size_t n3){
  _Complex double ***a;
  size_t n, m, k;

  a = malloc(n1*sizeof(_Complex double **));
  if ( a == NULL )
    sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
	      __LINE__);
    
  for ( n=0; n<n1; n++ ){
    a[n] = sep_complex_matrix(n2, n3);
    if ( a[n] == NULL )
      sep_error("%s at line %d: Couldn't allocate memory\n", __func__,
		__LINE__);
  }

  for ( n=0; n<n1; n++ )
    for ( m=0; m<n2; m++ )
      for ( k=0; k<n3; k++ )
	a[n][m][k] = 0.0 + 0.0*SEP_I;

  return a;
      
}

void sep_free_complex_matrix(_Complex double **a, size_t nr){
	
  for ( size_t i=0; i<nr; i++ ){
    free(a[i]);
  }
	
  free(a);
}


void sep_free_complex_tensor(_Complex double ***a, size_t n1, size_t n2){
  size_t n;
  
  for ( n=0; n<n1; n++ )
    sep_free_complex_matrix(a[n], n2);

  free(a);  

}

void sep_fprintf_complex_vector(_Complex double *a, const char *file, 
				char opt[2], size_t length){
  size_t i;
  
  FILE *fout = fopen(file, opt);
  if ( fout == NULL )
    sep_error("%s at %d: Couldn't open file.",__func__, __LINE__);

  for ( i=0; i<length; i++ )
    fprintf(fout, "%f ", creal(a[i]));
	
  fprintf(fout, "\n");
	
  for ( i=0; i<length; i++ )
    fprintf(fout, "%f ", cimag(a[i]));
	
  fprintf(fout, "\n");
	
  fclose(fout);
}

void sep_fprintf_complex_matrix(_Complex double **a, const char *file, 
				char opt[2], size_t nr, size_t nc){
  size_t i, j; 
	
  FILE *fout = fopen(file, opt);
  if ( fout == NULL )
    sep_error("%s at %d: Couldn't open file.",__func__, __LINE__);

  for ( i=0; i<nr; i++ ){
    for ( j=0; j<nc; j++ )
      fprintf(fout, "%f ", creal(a[i][j]));
    fprintf(fout, "\n");
  }
	
  for ( i=0; i<nr; i++ ){
    for ( j=0; j<nc; j++ )
      fprintf(fout, "%f ", cimag(a[i][j]));
    fprintf(fout, "\n");
  }
	
  fclose(fout);
}

#endif
