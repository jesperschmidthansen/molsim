/* 
* sepmisc.c - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/


#include "sepmisc.h"

void sep_error(char *str, ...){
  va_list ap;
  char *p, *sval, cval;
  int ival;
  double dval;
   
  printf("sep-error -> ");
   
  va_start(ap, str);
  for ( p=str; *p; p++ ){
    if ( *p != '%' ){
      putchar(*p);
      continue;
    }
    switch ( *++p ){
    case 'd': 
      ival = va_arg(ap, int);
      printf("%d", ival);
      break;
    case 'i': 
      ival = va_arg(ap, int);
      printf("%d", ival);
      break;
    case 'f':
      dval = va_arg(ap, double);
      printf("%f", dval);
      break;
    case 'c':
      cval = va_arg(ap, int);
      printf("%c", cval);
      break;
    case 's':
      for ( sval = va_arg(ap, char *); *sval; sval++ )	
	putchar(*sval);
      break;
    default:
      putchar(*p);
      break;
    }
  }

  printf("\n");
  printf("BAILING OUT\n");
  exit(EXIT_FAILURE);

}
		    

void sep_warning(char *str, ...){
  va_list ap;
  char *p, *sval, cval;
  int ival;
  double dval;
   
  printf("sep-warning -> ");
   
  va_start(ap, str);
  for ( p=str; *p; p++ ){
    if ( *p != '%' ){
      putchar(*p);
      continue;
    }
    switch ( *++p ){
    case 'd': 
      ival = va_arg(ap, int);
      printf("%d", ival);
      break;
    case 'i': 
      ival = va_arg(ap, int);
      printf("%d", ival);
      break;
    case 'f':
      dval = va_arg(ap, double);
      printf("%f", dval);
      break;
    case 'c':
      cval = va_arg(ap, int);
      printf("%c", cval);
      break;
    case 's':
      for ( sval = va_arg(ap, char *); *sval; sval++ )	
	putchar(*sval);
      break;
    default:
      putchar(*p);
      break;
    }
  }

  printf("\n");
  fflush(stdout);
}

double sep_lj(double r2, char opt){
  double retval=.0, rri, rri3;
  
  rri = 1.0/r2; rri3 = rri*rri*rri;
  
  switch (opt){
  case 'f':
    retval = 48. * rri3*(rri3 - 0.5)*rri;
    break;
  case 'u':
    retval = 4.0*rri3*(rri3-1.0); 
    break;
  }
  
  return retval;
}  

double sep_lj_shift(double r2, char opt){
  double retval=.0, rri, rri3;
  
  rri = 1.0/r2; rri3 = rri*rri*rri;
  
  switch (opt){
  case 'f':
    retval = 48.0 *rri3*(rri3 - 0.5)*rri;
    break;
  case 'u':
    retval = 4.0*rri3*(rri3 - 1.0) + SEP_LJCF2; 
    break;
  }
  
  return retval;
}  

double sep_wca(double r2, char opt){
  double retval=0.0, rri, rri3;
  
  rri = 1.0/r2; rri3 = rri*rri*rri;
  
  switch (opt){
  case 'f':
    retval = 48.0*rri3*(rri3 - 0.5)*rri;
    break;
  case 'u':
    retval = 4.0*rri3*(rri3 - 1.0) + 1.0; 
    break;
  }
  
  return retval;
}  


inline double sep_spring_x0(double r2, char opt){
  const double k = 500.0;
  double ret = 0.0;
  
  switch (opt){
  case 'f':
    ret = - k;
    break;
  case 'u':
    ret = 0.5*r2*k;
    break;
  }
  
  return ret;
}


int sep_reaction_1_order(seppart *ptr, const char *reaction, 
			 double crp, int npart){
  int n;
  unsigned int nreact;
  
  nreact = 0;
  for (n=0; n<npart; n++){
    if ((ptr[n].type==reaction[0]) && (sep_rand32() < crp)){
      nreact++;
      ptr[n].type = reaction[1];
    }
  }

  return nreact;
}


int sep_reaction_1_order_binomial(seppart *ptr, const char *reaction, 
			double crp, int npart){
  int molnumber;
  int nreact;
  
  int no_des_part = sep_count_type(ptr, reaction[0], npart);
  nreact = sep_binomial(crp,no_des_part);
  if(nreact==0) return nreact;
  int * v = sep_vector_int(nreact);
  for(int i=0;i<nreact;i++)
    {
      while(1)
	{
	  int already_flag=0;
	  molnumber = (int) (sep_rand32() * no_des_part);
	  for(int j=0;j<i;j++)
	    {
	      if(v[j]==molnumber) {already_flag=1;break;}
	    }
	  if(already_flag==0) break;
	}
      v[i]=molnumber;
    }
  qsort( v, nreact, sizeof(int), sep_compare_int_ascend);
  int moliter=0;
  int reaciter=0;
  for(int i=0;i<npart;i++)
    {
      if(ptr[i].type==reaction[0])
	{
	  if(moliter==v[reaciter])
	    {
	      ptr[i].type=reaction[1];
	      reaciter++;
	      if(reaciter==nreact) {free(v);return nreact;}
	    }
	  moliter++;
	}
    }
  return nreact;
}


unsigned int sep_reaction_2_order(seppart *ptr,  const char *rstr, 
				  double cr, double Pcr, double Q, 
				  sepsys *sys){ 
  int i1, i2, n, k;
  unsigned nreact = 0; 
  double r2, r[3], c_r2 = cr*cr, u, w, c; 
  static bool init = false;

  if ( init==false ){
    srand(time(NULL));
    init = true;
  }
    
  for (i1 = 0; i1 < sys->npart; i1++){ 
    n = 0; 
    while (1){ 
      i2 = ptr[i1].neighb[n]; 
      if (i2 == -1) break; 
      r2 = 0.0; 
      for (k=0; k<3; k++){ 
	r[k] = ptr[i1].x[k] - ptr[i2].x[k]; 
	sep_Wrap( r[k], sys->length[k] );
	r2 += r[k]*r[k];  
      }  
      if ( r2 < c_r2 && Pcr > sep_rand() ){  
	if ((ptr[i1].type == rstr[0]) &&  
 	    (ptr[i2].type == rstr[1])){  
 	  nreact++; 
 	  ptr[i1].type = rstr[2]; 
 	  ptr[i2].type = rstr[3]; 
          
          for ( k=0; k<3; k++ ){
            // Relativ vel: u = v_1 - v_2
            u = ptr[i1].v[k] - ptr[i2].v[k]; 
            // Reduced vel: w = (v_1 + v_2)/2
            w = (ptr[i1].v[k] + ptr[i2].v[k])*0.5; 
            // mu*u^2/2 + Q = mu*(c*u)^2/2 <=> 
            // c = (1 + 2*Q/mu*u^2)^(1/2), where mu is the 
            // reduced mass which is 0.5 here
            c = sqrt(1.0 + 4*Q/(u*u));
            // Post collision vel: 
            // v_1(post) = w + 0.5*U and  v_2(post) = w - 0.5*U 
            ptr[i1].v[k] = w + 0.5*c*u;
            ptr[i2].v[k] = w - 0.5*c*u;
          }
 	} 
 	else if ((ptr[i1].type == rstr[1]) &&  
 		 (ptr[i2].type == rstr[0])){ 
 	  nreact++; 
 	  ptr[i1].type = rstr[3]; 
 	  ptr[i2].type = rstr[2]; 
          
          for ( k=0; k<3; k++ ){
            // Relativ vel: u = v_1 - v_2
            u = ptr[i1].v[k] - ptr[i2].v[k]; 
            // Reduced vel: w = (v_1 + v_2)/2
            w = (ptr[i1].v[k] + ptr[i2].v[k])*0.5; 
            // mu*u^2/2 + Q = mu*(c*u)^2/2 <=> 
            // c = (1 + 2*Q/mu*u^2)^(1/2), where mu is the 
            // reduced mass which is 0.5 here
            c = sqrt(1.0 + 4*Q/(u*u));
            // Post collision vel: 
            // v_1(post) = w + 0.5*U and  v_2(post) = w - 0.5*U 
            ptr[i1].v[k] = w + 0.5*c*u;
            ptr[i2].v[k] = w - 0.5*c*u;
          }
 	} 
      } 
      n++;  
    }  
  }  
  
  return nreact; 
} 


void sep_save(seppart *ptr, size_t npart, const char *file){
  FILE *fout;

  fout = fopen(file, "wb");
  if ( fout == NULL )
    sep_error("%s at line %d: Couldn't open the file: %s\n", 
	      __func__, __LINE__, file);

  fwrite(ptr, sizeof(seppart), npart, fout);
  fclose(fout);

}


void sep_load(seppart *ptr, int nneighb, size_t npart, const char *file){
  FILE *fout;
  size_t n;

  for ( n=0; n<npart; n++ ) free(ptr[n].neighb);

  fout = fopen(file, "rb");
  if ( fout == NULL )
    sep_error("%s at line %d: Couldn't open the file: %s\n", 
	      __func__, __LINE__, file);

  if ( fread(ptr, sizeof(seppart), npart, fout) == 0 )
    sep_error("%s at line %d: Error reading file %s\n", 
	      __func__, __LINE__, file);
    
  fclose(fout);

  for ( n=0; n<npart; n++ ){
    ptr[n].neighb = malloc(nneighb*sizeof(int));
    if ( ptr[n].neighb == NULL ){
      sep_error("%s at line %d: Couldn't allocate memory\n",
		__func__, __LINE__);
    }
  }

}


void sep_relax_temp(seppart *ptr, char type, double Td, double tau,
		    sepsys *sys) {
  double fact;
  int n,k;

  double ekin = 0.0;  int ntype = 0;
  for ( n=0; n<sys->npart; n++ ){
    if ( ptr[n].type == type ){
      ntype ++;
      for ( k=0; k<3; k++ )
	ekin += ptr[n].m*sep_Sq(ptr[n].v[k]);
    }
  }

  ekin = 0.5*ekin;
  
  if ( ekin < DBL_EPSILON )
    sep_warning("%s at line %d: Zero kinetic energy - check your the types.",
		 __func__, __LINE__);
  
  double Ta = 2.0*ekin/(3*ntype);

  fact = sqrt(1.0+(sys->dt/tau)*(Td/Ta-1.0));

  for (n=0; n<sys->npart; n++){
    if ( ptr[n].type == type ){
      for (k=0; k<3; k++)
	ptr[n].v[k] *=  fact;
    }
  }

  sep_reset_momentum(ptr, type, sys);
  
}


void sep_reset_force(seppart *ptr, sepsys *sys){
  int n, k;

  for (n=0; n<sys->npart; n++)
    for (k=0; k<3; k++) ptr[n].f[k] = 0.0;

  sys->max_dist2 = 0.0;
}


void sep_reset_force_mol(sepsys *sys){
  unsigned int i, j;
  //  long int ii;

  if ( sys->molptr->flag_Fij == 0 )
    sep_error("%s: Tried to reset mol force, but flag is zero", __func__);

  for (  i=0; i<sys->molptr->num_mols; i++ ){
    for ( j=0; j<sys->molptr->num_mols; j++ ){
      for ( int k=0; k<3; k++ )  sys->molptr->Fij[i][j][k] = 0.0;
    }
  }

  /*
  //for ( i=0; i<(unsigned) sys->npart; i++ )
  for ( ii=0; ii< sys->npart; ii++ )
    for ( j=0; j<sys->molptr->num_mols; j ++ )
      for ( int k=0; k<3; k++ )  sys->molptr->Fiajb[ii][j][k] = 0.0;
  */
}


int sep_is_here(seppart *ptr, double min, 
		double max, int i, int dir){
  
  if ( (ptr[i].x[dir] > min) && (ptr[i].x[dir] < max) ) 
    return 1;
  else
    return 0;

}


void sep_get_xp(seppart*ptr, double dt, int npart, int ndim){
  int n, k;
  
  for (n=0; n<npart; n++)
    for (k=0; k<ndim; k++)
      ptr[n].xp[k] = ptr[n].x[k] - ptr[n].v[k]*dt;

}


int sep_nsubbox(double cf, double delta, double lbox){
  double cut;
  int nsubbox;

  cut = cf + delta;
  nsubbox = (int)(lbox/cut);
  
  return nsubbox;
}



double sep_box_length(double dens, int npart, int ndim){
  
  return  (pow((double)npart/dens, (double)1.0/ndim));
}

int sep_count_type(seppart *ptr,  char spec, int npart){
  int n, numb;

  numb = 0;
  for (n=0; n<npart; n++){
    if (ptr[n].type ==  spec)
      numb++;
  }

  return numb;
}


void sep_set_type(seppart *ptr,  char spec, int numb, int npart){
  int index, count;

  if ( sep_count_type(ptr, 'A', npart) < numb )
    sep_error("%s at line %d Need more 'n' particles\n",
	      __func__,__LINE__);
    
  srand(time(0));
  count=0;
  while (1){
    index = ceil(sep_rand()*npart)-1;
    if (ptr[index].type == spec)
      continue;
    else if (ptr[index].type == 'A'){
      ptr[index].type = spec;
      count++;
    }
    if (count == numb)
      break;
  }
      
}

void sep_set_x0(seppart *ptr, int npart){
  int n, k;
  
  for ( n=0; n<npart; n++ ){
    for ( k=0; k<3; k++ )
      ptr[n].x0[k] = ptr[n].x[k];
  }
  
}

void sep_set_xn(seppart *ptr, int npart){
  int n, k;
  
  for ( n=0; n<npart; n++ ){
    for ( k=0; k<3; k++ )
      ptr[n].xn[k] = ptr[n].x[k];
  }
  
}


void sep_save_xyz(seppart *ptr, const char *partnames, 
		  const char *file, char *mode, sepsys sys){
  FILE *fout;
  long int n, k, ntype, ntotal;
 
  ntype = strlen(partnames);
  
  fout = fopen(file, mode);
  if ( fout == NULL )
    sep_error("%s at line %d: I couldn't open file\n", __func__,
	      __LINE__);
  
  ntotal = 0;
  for ( n=0; n<ntype; n++ )
    ntotal += sep_count_type(ptr, partnames[n], sys.npart);
    
  fprintf(fout, "%lu\n%f %f %f\n", ntotal, 
	  sys.length[0], sys.length[1], sys.length[2]);

  for ( n=0; n<sys.npart; n++ ){
    for ( k=0; k<ntype; k++ ){
      if ( ptr[n].type == partnames[k] ){
	fprintf(fout, "%c %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n",  
		ptr[n].type , 
		ptr[n].x[0], ptr[n].x[1], ptr[n].x[2], 
		ptr[n].v[0], ptr[n].v[1], ptr[n].v[2], 
		ptr[n].m, ptr[n].z);
      }
    }
  }
  
  fclose(fout);
  
}


double sep_eval_mom(seppart *ptr, int npart){
  int n, k;
  double mom;

  mom = 0.0;
  for ( n=0; n<npart; n++ )
    for ( k=0; k<3; k++ )
      mom += ptr[n].v[k]*ptr[n].m;

  return mom/(npart*3);
    
}


double sep_eval_mom_type(seppart *ptr, char type, int dir, int npart){
  double mom;
  int n, ntype;

  mom = 0.0; ntype = 0;
  for (n=0; n<npart; n++){
    if (ptr[n].type == type){
      ntype ++;
      mom += ptr[n].v[dir]*ptr[n].m;
    }
  }
  
  return mom/ntype;
}

double sep_eval_momvec(seppart *ptr, int npart){
  int n, k;
  double momvec[3], mommag;

  for ( k=0; k<3; k++ )
      momvec[k] = 0.0;

  for ( n=0; n<npart; n++ )
    for ( k=0; k<3; k++ )
      momvec[k] += ptr[n].v[k]*ptr[n].m;

  mommag = 0.0;
  for ( k=0; k<3; k++ )
      mommag += momvec[k]*momvec[k];
  mommag = sqrt(mommag);

  return mommag;
    
}

void sep_set_density(seppart *ptr, char typef, char typet, 
                     double densf, double denst, int npart){
  unsigned int tt, count, i;
  
  tt = sep_count_type(ptr, typef, npart)*denst/densf;
    
  srand(42); count = 0;
  while ( 1 ){ 
    i = npart*sep_rand();
    if ( ptr[i].type == typef ){
      ptr[i].type = typet;
      count++;
    }
   
    if ( count == tt ) break;
  }
              
}


void sep_force_x0(seppart *ptr, char type, 
		  double (*fun)(double, char), sepsys *sys){
  double r2, r[3], f, ft, epot;
  int n, k;

  epot = .0; 
  for (n=0; n<sys->npart; n++){

    if ( ptr[n].type == type ){
      r2 = 0.0;
      for (k=0; k<3; k++){
	r[k] = ptr[n].x0[k] - ptr[n].x[k];
	sep_Wrap(r[k], sys->length[k] );
	r2 += r[k]*r[k];
      }
      ft = fun(r2, 'f');
      for (k = 0; k<3; k++){
	f = ft * r[k];
	ptr[n].f[k] -= f;
      }
      epot += fun(r2, 'u');
    }

  }

}




void sep_rand32_init(int *r250_index, int *r521_index, 
                     unsigned long *r250_buffer, 
                     unsigned long *r521_buffer){
  int i = SEP_R521_LEN;
  unsigned long mask1 = 1;
  unsigned long mask2 = 0xFFFFFFFF;
      
  while (i-- > SEP_R250_LEN) {
    r521_buffer[i] = rand();
  }
  while (i-- > 31) {
    r250_buffer[i] = rand();
    r521_buffer[i] = rand();
  }
  
  while (i-- > 0) {
    r250_buffer[i] = (rand() | mask1) & mask2;
    r521_buffer[i] = (rand() | mask1) & mask2;
    mask2 ^= mask1;
    mask1 >>= 1;
  }
  r250_buffer[0] = mask1;
  r521_buffer[0] = mask2;
  *r250_index = 0;
  *r521_index = 0;
}


double sep_rand32() {
#define R250_IA  (sizeof(unsigned long)*103)
#define R250_IB  (sizeof(unsigned long)*SEP_R250_LEN - R250_IA)
#define R521_IA  (sizeof(unsigned long)*168)
#define R521_IB  (sizeof(unsigned long)*SEP_R521_LEN - R521_IA)

  unsigned int count = 0;
  static int r250_index, r521_index; 
  static unsigned long r250_buffer[SEP_R250_LEN];
  static unsigned long r521_buffer[SEP_R521_LEN];
  int i1;
  int i2;
  unsigned char *b1 = (unsigned char*)r250_buffer;
  unsigned char *b2 = (unsigned char*)r521_buffer;
  unsigned long *tmp1, *tmp2;
  unsigned long r, s, a;
  int j1, j2;
  
  if ( count == 0 ){
    count = 1;
    sep_rand32_init(&r250_index, &r521_index, &r250_buffer[0], &r521_buffer[0]);
  }
  
  i1 = r250_index;
  i2 = r521_index;  
  b1 = (unsigned char*)r250_buffer;
  b2 = (unsigned char*)r521_buffer;
   
  j1 = i1 - R250_IB;
  if (j1 < 0)
    j1 = i1 + R250_IA;
  j2 = i2 - R521_IB;
  if (j2 < 0)
    j2 = i2 + R521_IA;
    
  tmp1 = (unsigned long *)(b1 + i1);
  r = (*(unsigned long *)(b1 + j1)) ^ (*tmp1);
  *tmp1 = r;
  tmp2 = (unsigned long *)(b2 + i2);
  s = (*(unsigned long *)(b2 + j2)) ^ (*tmp2);
  *tmp2 = s;
    
  i1 = (i1 != sizeof(unsigned long)*(SEP_R250_LEN-1)) ? 
    (i1 + sizeof(unsigned long)) : 0;
  r250_index = i1;
  i2 = (i2 != sizeof(unsigned long)*(SEP_R521_LEN-1)) ? (i2 + sizeof(unsigned long)) : 0;
  r521_index = i2;
  
  a = r^s;     
  
  return ((double)(a)/(UINT32_MAX))*2.0 - 1;
  
#undef R250_IA  
#undef R250_IB  
#undef R521_IA  
#undef R521_IB  

}

void sep_hs_update(sepatom *ptr, int ip, int jp, double tc, int flag,  
		   double lbox, size_t natoms, size_t ndim){
	
  size_t i, k;
  double r[3], v[3], b;
							      		  		 	
  // Updating the positions
  for ( i=0; i<natoms; i++){
    for ( k=0; k<ndim; k++ ){
      ptr[i].x[k] += ptr[i].v[k]*tc;
      sep_Periodic(ptr[i].x[k], lbox);                        
    }
  }                    
 
  // If collision occured find the new vel. of ip and jp
  if ( flag == 1 ) { 
    // Updating the velocities of ip and jp                  
    for ( k=0; k<ndim; k++ ){
      r[k] = ptr[ip].x[k] - ptr[jp].x[k];
      sep_Wrap(r[k], lbox);
      v[k] = ptr[ip].v[k] - ptr[jp].v[k];
    }
		
    b = sep_dot(r, v, ndim);
		
    for ( k=0; k<3; k++ ){
      ptr[ip].v[k] -= b*r[k];
      ptr[jp].v[k] += b*r[k];
    } 
  }
	
}

int sep_hscoll(int *ip, int *jp, double *tc,  seppart *ptr, int *list, 
	       sepsys *sys){
  double r[3], v[3], r2, v2, b, descr, tau;
  int j1, j2, m1, m1X, m1Y, m1Z, m2, m2X, m2Y, m2Z,
    k, offset, nsubbox3, nsubbox2, flag;
  static int iofX[] = {0,1,1,0,-1,0,1,1,0,-1,-1,-1,0,1};
  static int  iofY[] = {0,0,1,1,1,0,0,1,1,1,0,-1,-1,-1};
  static int  iofZ[] = {0,0,0,0,0,1,1,1,1,1,1,1,1,1};
  
  *tc = DBL_MAX;
  flag = 0;
  
  nsubbox2 = sys->nsubbox[1]*sys->nsubbox[0];
  nsubbox3 = sys->nsubbox[2]*nsubbox2;
   
	
  for (m1Z = 0; m1Z < sys->nsubbox[2]; m1Z++){
    for (m1Y = 0; m1Y < sys->nsubbox[1]; m1Y++) {
      for (m1X = 0; m1X < sys->nsubbox[0]; m1X++) {
	m1 = m1Z*nsubbox2 + m1Y*sys->nsubbox[0] + m1X;
	for (offset = 0; offset < 14; offset++) {
	 
	  m2X = m1X + iofX[offset];
	  if ( m2X == sys->nsubbox[0] ) {
	    m2X = 0;
	  } else if ( m2X == -1 ) {
	    m2X = sys->nsubbox[0]-1;
	  }
	  
	  m2Y = m1Y + iofY[offset];
	  if ( m2Y == sys->nsubbox[1] ) {
	    m2Y = 0;
	  } else if ( m2Y == -1 ) {
	    m2Y = sys->nsubbox[1]-1;
	  }
	  
	  m2Z = m1Z + iofZ[offset];
	  if ( m2Z == sys->nsubbox[2] ) {
	    m2Z = 0;
	  }
	  
	  j1 = list[m1];
	  m2 = m2Z*nsubbox2 + m2Y*sys->nsubbox[0] + m2X;
	  while ( j1 != -1 ) {
	    j2 = list[m2];
	    while ( j2 !=-1 ) {
	      if ( m1 != m2 || j2 < j1 ) {
		// asdf
		v2 = 0.0; r2 = 0.0;
		for ( k=0; k<3; k++ ){
		  r[k] = ptr[j1].x[k] - ptr[j2].x[k];
		  sep_Wrap(r[k], sys->length[k]);
		  v[k] = ptr[j1].v[k] - ptr[j2].v[k];
									
		  v2 += v[k]*v[k]; r2 += r[k]*r[k];
		}
	
		b = sep_dot(r, v, 3);
		descr = b*b - v2*(r2 - 1);
	
		if ( b < 0.0 && descr > 0.0 ){
		  tau =  (-b - sqrt(descr))/v2;
		  flag = 1;
		  if ( tau < *tc ){
		    *tc = tau;
		    *ip = j1;
		    *jp = j2;
		  }
		}
	      }
	      j2 = list[j2+nsubbox3];
	    }
	    j1 = list[j1+nsubbox3];
	  }
	}
      }
    }
  }
  
  return flag;
}


int sep_hs_coll(int *ip, int *jp, double *tc,  seppart *ptr, sepsys *sys){
  int *list, flag;
	
  list = sep_allocate_celllist(sys);
  sep_make_celllist(ptr, list, sys);
  flag = sep_hscoll(ip, jp, tc,  ptr, list, sys);

  free(list);

  return flag;
}



void sep_berendsen(sepatom *ptr, double Pd, double beta, sepret *ret, sepsys *sys){
  
  sep_pressure_tensor(ret, sys);

  double xi = 1 - beta*sys->dt*(Pd - ret->p);
  sys->length[2] *= xi;

  double scale = pow(xi, 1.0/3.0); 
  for ( int n=0; n<sys->npart; n++ ) ptr[n].x[2] *= scale;
  
  sys->volume = sys->length[0]*sys->length[1]*sys->length[2];
  
  // Need to update the sub boxes
  if ( sys->neighb_update != 0 ){
    sys->nsubbox[2] = sep_nsubbox(sys->cf, 0.0, sys->length[2]);
    if ( sys->nsubbox[2] < 3 ) 
      sep_warning("%s at %d: Number of subboxes in x direction are less than three",
		  __func__, __LINE__);
   
    sys->lsubbox[2] = sys->length[2]/sys->nsubbox[2];
  }

}
	


void sep_berendsen_iso(sepatom *ptr, double Pd, double beta, sepret *ret, sepsys *sys){


  sep_pressure_tensor(ret, sys);
  
  double xi = 1 - beta*sys->dt*(Pd - ret->p);

  for ( int k=0; k<3; k++ ) sys->length[k] *= xi;

  double scale = pow(xi, 1.0/3.0); 
  for ( int n=0; n<sys->npart; n++ )
    for ( int k=0; k<3; k++ ) ptr[n].x[k] *= scale;
  
  sys->volume = sys->length[0]*sys->length[1]*sys->length[2];
  
  // Need to update the sub boxes
  if ( sys->neighb_update != 0 ){
    for ( int k=0; k<3; k++ ){
      sys->nsubbox[k] = sep_nsubbox(sys->cf, 0.0, sys->length[k]);
      if ( sys->nsubbox[k] < 3 ) 
	sep_warning("%s at %d: Number of subboxes in x direction are less than three",
		    __func__, __LINE__);
   
      sys->lsubbox[k] = sys->length[k]/sys->nsubbox[k];
    }
  }
  
}	


void sep_berendsen_mol(sepatom *ptr, sepmol *mol, double Pd, 
		       double beta, sepret *ret, sepsys *sys){


  sep_mol_cm(ptr, mol, sys);
  double *dr0 = sep_vector(sys->npart);
  for ( int n=0; n<sys->npart; n++ ){
    int i = ptr[n].molindex;
    if ( i != -1 ) dr0[n] = ptr[n].x[2] - mol[i].x[2];
    
  }
  
  sep_mol_pressure_tensor(ptr,  mol, ret, sys);

  double xi = 1 - beta*sys->dt*(Pd - ret->p_mol);
  sys->length[2] *= xi;
  sys->volume = sys->length[0]*sys->length[1]*sys->length[2];

  double scale = pow(xi, 1.0/3.0); 
  for ( unsigned n=0; n<sys->molptr->num_mols; n++ ) 
    mol[n].x[2] *= scale;

  for ( int n=0; n<sys->npart; n++ ){
    int i = ptr[n].molindex;
    if ( i != -1 ) {
      ptr[n].x[2] = mol[i].x[2] + dr0[n];
      sep_Periodic( ptr[n].x[2], sys->length[2] );
    }
  }

  // In case we have neighbour-list we need to update the sub boxes
  if ( sys->neighb_update != 0 ){
    sys->nsubbox[2] = sep_nsubbox(sys->cf, 0.0, sys->length[2]);
    if ( sys->nsubbox[2] < 3 ) 
      sep_warning("%s at %d: Number of subboxes in x direction are less than three",
		  __func__, __LINE__);
    sys->lsubbox[2] = sys->length[2]/sys->nsubbox[2];
  }


  free(dr0);
  
}	



void sep_compress_box(sepatom *ptr, double rhoD, double xi, sepsys *sys){

  const double density = sys->npart/sys->volume; 

  if ( fabs(density - rhoD) < 1e-6 ) return; // ACHTUNG
  
  if ( density > rhoD ) xi = 1.0/xi;
    
  for ( int k=0; k<3; k++){
    sys->length[k] *= xi;
    if  ( sys->length[k] < sys->cf*2.0 )
      sep_warning("%s at %d: Box length too small compared to the maximum cut-off", 
		  __func__, __LINE__);
  }

  for ( int n=0; n<sys->npart; n++ )
    for ( int k=0; k<3; k++ ) ptr[n].x[k] *= xi;
    
  // Need to update the sub boxes
  if ( sys->neighb_update != 0 ){

    for ( int k=0; k<3; k++ ){
      sys->nsubbox[k] = sep_nsubbox(sys->cf, sys->skin, sys->length[k]);
      if ( sys->nsubbox[k] < 3 ) 
	sep_warning("%s at %d: Number of subboxes in x direction are less than three", 
		    __func__, __LINE__);
      sys->lsubbox[k] = sys->length[k]/sys->nsubbox[k];
    }
  }
    
  sys->volume = sys->length[0]*sys->length[1]*sys->length[2];

}


void sep_compress_box_dir(sepatom *ptr, double rhoD, double xi, 
			  int dir, sepsys *sys){

  if ( sys->npart/sys->volume < rhoD ){ 
    
    sys->length[dir] *= xi;
    if  ( sys->length[dir] < sys->cf*2.0 )
      sep_warning("%s at %d: Box length too small compared to the maximum cut-off", 
		   __func__, __LINE__);
			   
    double scale = pow(xi, 1.0/3.0); 
    for ( int n=0; n<sys->npart; n++ ) ptr[n].x[dir] *= scale;
    
    // Need to update the sub boxes
    if ( sys->neighb_update != 0 ){
      sys->nsubbox[dir] = sep_nsubbox(sys->cf, 0.0, sys->length[dir]);
      if ( sys->nsubbox[dir] < 3 ) 
	sep_warning("%s at %d: Number of subboxes in x direction are less than three", 
		    __func__, __LINE__);
    
      sys->lsubbox[dir] = sys->length[dir]/sys->nsubbox[dir];
    }

    sys->volume = sys->length[0]*sys->length[1]*sys->length[2];
  }

}

void sep_compress_box_dir_length(sepatom *ptr, double length, double xi, 
				 int dir, sepsys *sys){

  if ( sys->length[dir] > length ){ 
    
    sys->length[dir] *= xi;
    if  ( sys->length[dir] < sys->cf*2.0 )
      sep_warning("%s at %d: Box length too small compared to the maximum cut-off", 
		   __func__, __LINE__);
			   
    double scale = pow(xi, 1.0/3.0); 
    for ( int n=0; n<sys->npart; n++ ) ptr[n].x[dir] *= scale;
    
    // Need to update the sub boxes
    if ( sys->neighb_update != 0 ){
      sys->nsubbox[dir] = sep_nsubbox(sys->cf, 0.0, sys->length[dir]);
      if ( sys->nsubbox[dir] < 3 ) 
	sep_warning("%s at %d: Number of subboxes in x direction are less than three", 
		    __func__, __LINE__);
    
      sys->lsubbox[dir] = sys->length[dir]/sys->nsubbox[dir];
    }

    sys->volume = sys->length[0]*sys->length[1]*sys->length[2];
  }

}



void sep_set_charge(seppart *ptr, char type, double z, sepsys sys){

  for ( int n=0; n<sys.npart; n++ )
    if ( ptr[n].type == type ) ptr[n].z = z;

}

void sep_set_mass(seppart *ptr, char type, double m, sepsys sys){

  for ( int n=0; n<sys.npart; n++ )
    if ( ptr[n].type == type ) ptr[n].m = m;

}

double sep_dist_ij(double *r, seppart *ptr, int i, int j, sepsys *sys){


  double r2 = 0.0;
  for ( int k=0; k<3; k++ ){
    r[k]  = ptr[i].x[k]-ptr[j].x[k];
    sep_Wrap( r[k], sys->length[k] );
    r2   += r[k]*r[k];
  }

  return sqrt(r2);
}



void sep_set_omp(unsigned nthreads, sepsys *sys){
	
  sys->omp_flag = true;
  sys->nthreads = nthreads;

  omp_set_num_threads(nthreads);

}

void sep_set_skin(sepsys *sys, double value){
	
  sys->skin = value;
	
}

double sep_randn(void){
  static unsigned int flag = 0;
  static double r1=0.0;
  static double r2=0.0;

  double w, x1, x2;

  if ( flag == 1 ){
    flag = 0;
    return r2;
  }
  else {
    do {
      x1 = 2.0 * sep_rand() - 1.0;
      x2 = 2.0 * sep_rand() - 1.0;
      
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 || w == 0.0 );
  
    w = sqrt( (-2.0 * log( w ) ) / w );
    r1 = x1 * w;
    r2 = x2 * w;

    flag = 1;

    return r1;
  }

}


void sep_set_ldiff(sepatom *ptr, char type, double ldiff, sepsys sys){

  for ( int n=0; n<sys.npart; n++ ){
    if ( ptr[n].type == type )
      ptr[n].ldiff = ldiff;

  }

}


void sep_reset_momentum(seppart *ptr, const char type, sepsys *sys){
  double sum_mom[3]={0.0};

  double mass = 0.0; int ntype = 0;
  for ( int n=0; n<sys->npart; n++ ){
    if ( ptr[n].type == type ){
      ntype++;
      for ( int k=0; k<3; k++ ){
	sum_mom[k]  += ptr[n].v[k]*ptr[n].m;
      }
      mass += ptr[n].m;
    }
  }
  
  for ( int n=0; n<sys->npart; n++ )
    if ( ptr[n].type == type ) 
      for ( int k=0; k<3; k++ )
	ptr[n].v[k] -=  sum_mom[k]/mass;

}

	
	
void sep_set_ndof(size_t ndof, sepsys *sys){

  sys->ndof = ndof;

}

void sep_eval_xtrue(seppart *ptr, sepsys *sys){

  for ( int n=0; n<sys->npart; n++ )
    for ( int k=0; k<3; k++ )
      ptr[n].xtrue[k] = sys->length[k]*ptr[n].crossings[k] + ptr[n].x[k];

}



double sep_ran0(long *idum){

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
  
  long k;
  double ans;
  
  *idum ^= MASK;
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  ans=AM*(*idum);
  *idum ^= MASK;
  
  return ans;

#undef IA 
#undef IM 
#undef AM 
#undef IQ 
#undef IR 
#undef MASK 

}



double sep_ran3(long *idum){

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

  static int inext,inextp;
  static long ma[56]; 
  static int iff=0; 
  long mj,mk;
  int i,ii,k;
  if (*idum < 0 || iff == 0) { 
    iff=1;
    mj=labs(MSEED-labs(*idum)); 
    mj %= MBIG; 
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) { 
      ii=(21*i) % 55; 
      ma[ii]=mk; 
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++) 
      for (i=1;i<=55;i++) { 
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0; 
    inextp=31; 
    *idum=1;
  }
  
  if (++inext == 56) inext=1; 
  if (++inextp == 56) inextp=1; 
  mj=ma[inext]-ma[inextp]; 
  if (mj < MZ) mj += MBIG; 
  ma[inext]=mj; 
  return mj*FAC; 

#undef MBIG 
#undef MSEED 
#undef MZ 
#undef FAC 

}
