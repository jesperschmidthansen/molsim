
/* 
* seputil.c - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/

#include "seputil.h"

void sep_solvelineq(double *x, double **A, int n){
  double *vMaxI, v, vMax, tmp;
  int *ptrMax; 
  int i, j, k, m;

  vMaxI =sep_vector(n);
  ptrMax = sep_vector_int(n);
  for ( i=0; i<n; i++ )
    ptrMax[i] = 0;

  for (i = 0; i < n; i ++) {
    vMax = 0.;
    for (j = 0; j < n; j ++) {
      if ((v = fabs(A[i][j])) > vMax) vMax = v;
    }
    vMaxI[i] = 1. / vMax;
  }
  for (m = 0; m < n; m ++) {
    vMax = 0.;
    for (i = m; i < n; i ++) {
      for (k = 0; k < m; k ++) A[i][m] -= A[i][k] * A[k][m];
      if ((v = fabs (A[i][m]) * vMaxI[i]) > vMax) {
        vMax = v;
        ptrMax[m] = i;
      }
    }
    if (m != ptrMax[m]) {
      for (k = 0; k < n ; k ++){
	tmp = A[m][k];
	A[m][k] =  A[ptrMax[m]][k];
	A[ptrMax[m]][k] = tmp;
      }
      vMaxI[ptrMax[m]] = vMaxI[m];
    }
    for (j = m + 1; j < n; j ++) {
      for (k = 0; k < m; k ++) 
	A[m][j] -= A[m][k] * A[k][j];
      A[m][j] /= A[m][m];
    }
  }
  for (i = 0; i < n; i ++) {
    tmp = x[ptrMax[i]];
    x[ptrMax[i]] = x[i];
    x[i] = tmp;
    for (j = 0; j < i; j ++) 
      x[i] -= A[i][j] * x[j];
    x[i] /= A[i][i];
  }

  for (i = n - 2; i >= 0; i --) {
    for (j = i + 1; j < n; j ++) 
      x[i] -= A[i][j] * x[j];
  }

  free(ptrMax); free(vMaxI);
}


void sep_solvelineq1(double *ret, double *x, double **A, int n){
  double *vMaxI, v, vMax, tmp, *tmpx, **tmpA;
  int *ptrMax; 
  int i, j, k, m;

  // Copy to save one return
  tmpx = sep_vector(n);
  tmpA = sep_matrix(n, n);
  for ( i=0; i<n; i++ ){
    tmpx[i] = x[i];
    for ( j=0; j<n; j++ )
      tmpA[i][j] = A[i][j];
  }

  vMaxI =sep_vector(n);
  ptrMax = sep_vector_int(n);
  for ( i=0; i<n; i++ )
    ptrMax[i] = 0;

  for (i = 0; i < n; i ++) {
    vMax = 0.;
    for (j = 0; j < n; j ++) {
      if ((v = fabs(A[i][j])) > vMax) vMax = v;
    }
    vMaxI[i] = 1. / vMax;
  }
  for (m = 0; m < n; m ++) {
    vMax = 0.;
    for (i = m; i < n; i ++) {
      for (k = 0; k < m; k ++) A[i][m] -= A[i][k] * A[k][m];
      if ((v = fabs (A[i][m]) * vMaxI[i]) > vMax) {
        vMax = v;
        ptrMax[m] = i;
      }
    }
    if (m != ptrMax[m]) {
      for (k = 0; k < n ; k ++){
	tmp = A[m][k];
	A[m][k] =  A[ptrMax[m]][k];
	A[ptrMax[m]][k] = tmp;
      }
      vMaxI[ptrMax[m]] = vMaxI[m];
    }
    for (j = m + 1; j < n; j ++) {
      for (k = 0; k < m; k ++) 
	A[m][j] -= A[m][k] * A[k][j];
      A[m][j] /= A[m][m];
    }
  }
  for (i = 0; i < n; i ++) {
    tmp = x[ptrMax[i]];
    x[ptrMax[i]] = x[i];
    x[i] = tmp;
    for (j = 0; j < i; j ++) 
      x[i] -= A[i][j] * x[j];
    x[i] /= A[i][i];
  }

  for (i = n - 2; i >= 0; i --) {
    for (j = i + 1; j < n; j ++) 
      x[i] -= A[i][j] * x[j];
  }


  for ( i=0; i<n; i++ ){
    ret[i] = x[i];
    x[i] = tmpx[i];
    for ( j=0; j<n; j++ )
      A[i][j] = tmpA[i][j];
  }

  free(ptrMax); free(vMaxI); free(tmpx);
  sep_free_matrix(tmpA, n);

}


void sep_linfit(double *param, double *x, double *data,  
		size_t ndata,  size_t nparam, double (*fun)(double, int)){
  double **a, *b, yi;
  size_t k, j, i;

  a = sep_matrix(ndata, ndata);
  b = sep_vector(nparam);

  for ( k=0; k<nparam; k++ ){
    b[k] = 0.0;
    for ( j=0; j<nparam; j++ ){
      a[k][j] = 0.0;
      for ( i=0; i<ndata; i++ ){
	yi = fun(x[i], k);
	a[k][j] += fun(x[i], j)*yi;
	b[k] += data[i]*yi;
      }
    }
  }

  sep_solvelineq(b, a, nparam);

  for ( k =0; k<nparam; k++ ) param[k] = b[k]/nparam; 

  sep_free_matrix(a, ndata);
  free(b);

}


void sep_linfit_list(double *param, double *x, 
		     double *data, size_t ndata,  
		     unsigned int *paramlist, size_t nparam,
		     double (*fun)(double, int)){
  double **a, *b, yi;
  size_t k, j, i;

  a = sep_matrix(ndata, ndata);
  b = sep_vector(nparam);

  for ( k=0; k<nparam; k++ ){
    b[k] = 0.0;
    for ( j=0; j<nparam; j++ ){
      a[k][j] = 0.0;
      for ( i=0; i<ndata; i++ ){
				yi = fun(x[i], paramlist[k]);
				a[k][j] += fun(x[i], paramlist[j])*yi;
				b[k] += data[i]*yi;
      }
    }
  }

  sep_solvelineq(b, a, nparam);

  for ( k=0; k<nparam; k++ ) param[k] = b[k]/nparam;

  sep_free_matrix(a, ndata);
  free(b);

}


double sep_vector_mean(double *vec, int length){
  int n; 
  double sum;

  sum = 0.;
  for (n=0; n<length; n++)
    sum += vec[n];

  return (sum/length);
}


double sep_vector_std(double *vec, int length){
  int n;
  double sum, mean, var;

  sum = 0.;  var = 0.;
  for (n=0; n<length; n++)
    sum += vec[n];

  mean = sum/length;

  for (n=0; n<length; n++)
    var += (vec[n] - mean)*(vec[n]-mean);

  return sqrt(var/(length-1));
}



double sep_vector_kurtosis(double *vec, int length){
  int n; 
  double v, v4, sumv4, sumv2, mean, tl, nv;
  
  mean = sep_vector_mean(vec, length);

  sumv4 = sumv2 = 0.0;
  for (n=0; n<length; n++){
    v = vec[n] - mean;                            
    v4 = pow(v, 4.0);
    sumv4 += v4;
    sumv2 += pow(v, 2.0);
  }

  tl = sumv4/length;
  nv = pow(sumv2/(length - 1.0), 2.0);

  return (tl/nv - 3.0); 
}

double sep_vector_skewness(double *vec, int length){
  int n;  
  double v, v3, sumv3, sumv2, mean, tl, nv; 
  
  mean = sep_vector_mean(vec, length);

  sumv3 = sumv2 = 0.0;
  for (n=0; n<length; n++){
    v = vec[n] - mean;                            
    v3 = pow(v, 3.0);
    sumv3 += v3;
  }
  
  tl = sumv3/length;
  nv = pow(sep_vector_std(vec, length), 3.0);

  return (tl/nv); 
}


double sep_vector_max(double *vec, int length){
  int n;
  double max;

  max = vec[0];
  for (n=1; n<length; n++){
    if (vec[n] > max) max = vec[n];
  }

  return max;
}


double sep_vector_min(double *vec, int length){
  int n;
  double min;

  min = vec[0];
  for (n=1; n<length; n++){
    if (vec[n] < min) min = vec[n];
  }

  return min;
} 


void sep_autocorr(double *out, double *in, size_t lvec){
  size_t i, j;
	
  for ( i=0; i<lvec; i++ ){
    out[i] = 0.0;	
    for ( j=0; j<lvec-i; j++ )
      out[i] += in[j]*in[j+i];	
  }
	
  for ( i=0; i<lvec; i++ )
    out[i] /= (lvec-i);		
}


double sep_lincorr(double *x, double *y, size_t ndata){
  double sumx, sumy, meanx, meany, sumdxy, sumdx2, sumdy2, dx, dy;
  size_t i;

  sumx = 0.0;
  sumy = 0.0;
  for ( i=0; i<ndata; i++ ){
    sumx += x[i];
    sumy += y[i];
  }

  meanx = sumx/ndata;
  meany = sumy/ndata;

  sumdxy = 0.0; 
  sumdx2 = 0.0; 
  sumdy2 = 0.0;
  for ( i=0; i<ndata; i++ ){
    dx = x[i]- meanx;
    dy = y[i] - meany;
    sumdx2 += dx*dx;
    sumdy2 += dy*dy;
    sumdxy += dx*dy;
  }

  return sumdxy/sqrt(sumdx2*sumdy2);

}


double sep_trapz(double *f, double x_0, double x_1, int length){
  double s, h;
  int n;

  h = (x_1-x_0)/length;

  s = 0.5*f[0];
  for (n=1; n<length-1; n++)
    s += f[n];
  s += 0.5*f[length-1];

  return s*h;
}



 int sep_sgn(double x){

  // Asuming x == 0.0 never happens!
  if ( x < 0.0 ){
    return -1;
  }
  else {
    return 1;  
  }
}


 int sep_heaviside(double x){

  // Asuming x==0.0 never happens!
  if ( x < 0.0 ){
    return 0;
  }
  else {
    return 1;
  }

}

 double sep_dot(double *a, double *b, int length){
  double retval;
  int n;

  retval = 0.0;
  for ( n=0; n<length; n++ )
    retval += a[n]*b[n];
  
  return retval;

}

 double sep_vector_norm(double *a, int length){
  int n;
  double norm = 0.0;
  
  for ( n=0; n<length; n++ )  norm += a[n]*a[n];
    
  return sqrt(norm);

}

double sep_trace(double **A, size_t nrow){
  double trace = 0.0;
  size_t n;
  
  for ( n=0; n<nrow; n++ )
    trace += A[n][n]; 	
   	
  return trace;      
}

void sep_vector_multiply(double *a, double value, size_t lvec){
  size_t n; 
  
  for ( n=0; n<lvec; n++ ) a[n] *= value;

}

void sep_vector_add(double *a, double value, size_t lvec){
  size_t n; 
  
  for ( n=0; n<lvec; n++ ) a[n] += value;

}

void sep_matrix_multiply(double **A, double value, 
                                size_t nrow, size_t ncol){
  size_t n, m; 
  
  for ( n=0; n<nrow; n++ )
    for ( m=0; m<ncol; m++ ) A[n][m] *= value;

}

void sep_matrix_add(double **A, double value, 
                           size_t nrow, size_t ncol ){
  size_t n, m; 
  
  for ( n=0; n<nrow; n++ ) 
    for ( m=0; m<ncol; m++ ) A[n][m] += value;

}


void sep_cross3(double *ret, double *a, double *b){
  
  ret[0] = a[1]*b[2] - a[2]*b[1];
  ret[1] = a[2]*b[0] - a[0]*b[2];
  ret[2] = a[0]*b[1] - a[1]*b[0];

}

void sep_linspace(double *vec, double a, double b, size_t lvec){
  double delta = (b-a)/(lvec-1);			
	
  vec[0] = a;
  for ( size_t i=1; i<lvec; i++ )
    vec[i] = vec[i-1] + delta;
	
}


void sep_decompose_tensor(double *p, double **Pos, double **Pa, double **P, 
			  size_t nrc){
  size_t i, j;
  *p = sep_trace(P, nrc)/3.0;
	
  for ( i=0; i<nrc; i++ ){
    for ( j=0; j<nrc; j++ ){
      Pa[i][j] = 0.5*(P[i][j] - P[j][i]);
      Pos[i][j] = 0.5*(P[i][j] + P[j][i]);
      if ( j==i )
	Pos[i][j] -= *p;	
    }
  }
	
}
		
void sep_mv_multi(double *out, double **A, double *x, int nr, int nc){

  for ( int n = 0; n<nr; n++ ){
    out[n] = 0.0;
    for ( int m = 0; m<nc; m++ )
      out[n] += A[n][m]*x[m];
  }

}


unsigned sep_maxeigvalue_pow(double *q, double **A, int nrc, double eps0, 
			     const unsigned countermax){
  double x1[3]={1.0}, x2[3], eps, maxx, m1, m2, m3;
  int n;
  unsigned counter = 0;

  while ( 1 ){
    
    sep_mv_multi(x2, A, x1, nrc, nrc);
    
    m1 = 0.0;
    m2 = 0.0;
    m3 = 0.0;
    for ( n=0; n<nrc; n++ ){
      m1 += x1[n]*x1[n];
      m2 += x1[n]*x2[n];
      m3 += x2[n]*x2[n];
    }

    *q = m2/m1;
   
    eps = sqrt(m3/m1 - (*q)*(*q));
    counter ++;
    if ( eps < eps0 || counter == countermax ) break;

    maxx =0.0;
    for ( n=0; n<nrc; n++ ){
      if ( fabs(x2[n]) > maxx ) maxx = x2[n];
    }

    for ( n=0; n<nrc; n++ ) x1[n] = x2[n]/maxx;
  }
   
  return counter;
}
		
		
double sep_vector_sum(double *a, size_t length){
  size_t i;
  double sum = 0;

  for ( i=0; i<length; i++ ) sum += a[i];

  return sum;
}

		
double sep_det3(double **A){
  
  double a = A[0][0]*A[1][1]*A[2][2];
  double b = A[0][0]*A[1][2]*A[2][1];
  double c = A[0][1]*A[1][0]*A[2][2];
  double d = A[0][1]*A[1][2]*A[2][0];
  double e = A[0][2]*A[1][0]*A[2][1];
  double f = A[0][2]*A[1][1]*A[2][0];

  return a - b - c + d + e - f;

}
		
int sep_binomial(double p, int n)
/* using inverese cdf method, crushes if n*p is large
 * (then it's better to approximate with a gaussian)
 */

{
  int ret;
  double q;
  double ratio;
  long double f;
  double cdf;
  
  q = 1. - p;
  ratio = p / q;

  cdf = sep_rand32();
  f = powl(q,n);

  for(ret = 0; ret <= n; ++ret)
    {
      if(cdf < f)
	return ret;
      cdf -= f;
      f *= ratio * (n-ret) / (ret + 1);
    }
  return ret;
}

int sep_compare_int_ascend( const void *a, const void *b)
{
     int int_a = * ( (int*) a );
     int int_b = * ( (int*) b );

     if ( int_a < int_b ) return -1;
     else if ( int_a > int_b ) return 1;
     else return 0;
}


void sep_eig_real_symmetric(double *eig, double** A){
  double phi, q, p, r, p1, p2, **B;

  B = sep_matrix(3,3);
  
  p1 = sep_Sq(A[0][1]) + sep_Sq(A[0][2]) + sep_Sq(A[1][2]);
  if ( p1 == 0.0 ){     // A is diagonal.
    eig[0] = A[0][0];
    eig[1] = A[1][1];
    eig[2] = A[2][2];
  }
  else {
    q = sep_trace(A, 3)/3.0;
    p2 = sep_Sq(A[0][0] - q) + sep_Sq(A[1][1] - q) + 
      sep_Sq(A[2][2] - q) + 2*p1;
    p = sqrt(p2/6.);
   
    for ( int i=0; i<3; i++ ){
      for ( int j=0; j<3; j++ ){
	if ( j==i )
	  B[i][j] = 1./p*(A[i][j] - q);
	else 
	  B[i][j] = 1./p*A[i][j];
      }
    }
    r = sep_det3(B)/2.0;

    if (r <= -1.0 ) 
      phi = SEP_PI/3.0;
    else if (r >= 1.0 )
      phi = 0.0;
    else
      phi = acos(r)/3.0;
	
    eig[0] = q + 2.0*p*cos(phi);
    eig[1] = q + 2.0*p*cos(phi + (2*SEP_PI/3.0));
    eig[2] = 3.0*q - eig[0] - eig[1];
    
  }

  sep_free_matrix(B, 3);
}



void sep_corr_complex_arrays(_Complex double **A, _Complex double **B,
                             _Complex double **C, size_t nrows, size_t ncol){

  size_t k, n, nn;

#pragma omp parallel for                        \
  private(n, nn)                                
  for ( k=0; k<ncol; k++ )
    for ( n=0; n<nrows; n++ )
      for ( nn=0; nn<nrows-n; nn++ ) A[n][k] += B[nn][k]*C[n+nn][k];
  

}
