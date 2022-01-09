#include "septools.h"

 
void allocate_memory(data *dptr){
  const unsigned numb = dptr->npart_type;

  dptr->m = sep_vector(numb);
  dptr->z = sep_vector(numb);

  dptr->x = sep_matrix(numb, 3);
  dptr->v = sep_matrix(numb, 3);
  dptr->s = sep_matrix(numb, 3);
  dptr->w = sep_matrix(numb, 3);
  dptr->inertia = sep_matrix(numb, 6);

}

void free_memory(data *dptr){
  const unsigned numb = dptr->npart_type;

  free(dptr->m); 
  free(dptr->z);

  sep_free_matrix(dptr->x, numb);
  sep_free_matrix(dptr->v, numb);
  sep_free_matrix(dptr->s, numb);
  sep_free_matrix(dptr->w, numb);
  sep_free_matrix(dptr->inertia, numb);

}


int file_exists(const char file[]){

  FILE *fp = fopen(file, "r");
  if( fp ) {
    fclose(fp);
    return 1;
  } 
  else 
    return 0;
}


void get_num_type(data *dptr, const char type, const char file[]){

  double tnow, lx, ly, lz;
  char line[2048];

  dptr->npart = 0; 
  dptr->npart_type = 0;
  dptr->type = type;

  FILE *fin = fopen(file, "r");
  if ( fin == NULL ){
    sep_error("Couldn't read input file- BAILING OUT");
  }

  fscanf(fin, "%lf %lf %lf %lf\n", &tnow, &lx, &ly, &lz);
  
  while ( fgets(line, 2048, fin) != NULL ) {
    (dptr->npart)++;
    if ( line[0] == type ) (dptr->npart_type)++;
  }
  
  fclose(fin);

  printf("Found %d molecules of type %c (total number of molecules where %d)\n",
	 dptr->npart_type, type, dptr->npart);
  fflush(stdout);

}

void get_time(data *dptr, const char file[]){

  double tnow, lx, ly, lz;

  FILE *fin = fopen(file, "r");
  if ( fin == NULL ){
    sep_error("Couldn't read input file- BAILING OUT");
  }

  fscanf(fin, "%lf %lf %lf %lf\n", &tnow, &lx, &ly, &lz);

  fclose(fin);
  dptr->tnow = tnow;
}


void read_entry(data *dptr, const char file[]){
  double tnow, lx, ly, lz;
  char type;
  int nuau;
  double m, x[3], v[3], s[3], w[3], inertia[6]; // xx,yy,zz,xy,xz,yz component
  
  FILE *fin = fopen(file, "r");
  if ( fin == NULL )
    sep_error("Couldn't read input file- BAILING OUT");
  
  fscanf(fin, "%lf %lf %lf %lf\n", &tnow, &lx, &ly, &lz);
  
  int j = 0;
  for ( unsigned i=0; i<dptr->npart; i++){
    int nread = fscanf(fin, "%c %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
		       &type, &nuau, &m, 
		       &x[0], &x[1], &x[2], 
		       &v[0], &v[1], &v[2],
		       &s[0], &s[1], &s[2],
		       &w[0], &w[1], &w[2],
		       &inertia[0], &inertia[1], &inertia[2],
		       &inertia[3], &inertia[4], &inertia[5]);
    if ( nread != 21 )
      sep_error("Read failure - BAILING OUT");

    if ( type == dptr->type ){
      dptr->m[j] = m;
      for ( int k=0; k<3; k++ ) {
	dptr->x[j][k] = x[k];
	dptr->v[j][k] = x[k];
	dptr->s[j][k] = x[k];
	dptr->w[j][k] = x[k];
      }
      for ( int k=0; k<6; k++ ) {
	dptr->inertia[j][k] = inertia[k];
      }
      j++;
    }

  }
  
  fclose(fin);
} 
