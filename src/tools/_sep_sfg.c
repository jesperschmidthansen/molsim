
/*************************************************
 * 
 * sep_sfg.c
 *
 * Generates xyz and top files from individual config files
 *
 * Copyright, 2021,  Jesper Schmidt Hansen
 *
 *******************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define VERSION "0.1"

#define RAND() ( rand()/(RAND_MAX+1.0) )


/************************** STRUCTS ****************************/

// Represents a molecule defined via  
// the xyz and top files
typedef struct {

  unsigned nmol_type, nuau;

  float *rel_cm;
  char *types;
  float *mass;
  float *z;

  unsigned *bonds, nbonds;
  unsigned *angles, nangles;
  unsigned *dihedrals, ndihedrals;

  float cm[3];

} molecule;



// Lattice on which the molecules are placed
typedef struct {
  float *x, *y, *z;
  float lbox;
} lattice;



/**************** DECLARATIONS/PROTOTYPES *******************/

// Aux. function
void merror(const char *str, const char *fun, int line);
FILE *setfilepointer(FILE *fptr, const char *section);
void pipe( molecule *mol, int ntypes );
float velocity( void );
int findnumtypes(molecule *mol, int ntypes_mol);
void printpre(int nmol_tot, int ntypes );
void  printinsect(int n, int ntypes);
float wrap(float r, float lbox);

// Memory management 
molecule *allmem(unsigned nmoltypes);
lattice setcubic(unsigned nmol_tot, float dens);
void freemem(lattice *lptr, molecule *mptr, unsigned nmols, unsigned *rlist);

// Postion information 
unsigned readxyz(molecule *mol, const char *file, unsigned moli, unsigned nmol);
void writexyz(molecule *mol, unsigned moli,  lattice lat, unsigned  lati, 
	      float lbox,  int seed, const char *str);

void writexyzheader(unsigned natoms, float lbox, const char *str);
unsigned *setmolpos( unsigned ntypes, unsigned nmol, unsigned *ntypes_per_mol, 
		     unsigned seed );

// Bond information
void readtop_bonds(molecule *mol, const char *fstr, unsigned molt);
void writebondheader(const char *fstr);
void writebondssect(molecule *mol, unsigned molt, unsigned moli, const char *fstr);


// Angle information
void readtop_angles(molecule *mol, const char *fstr, unsigned molt);
void writeangleheader(const char *fstr);
void writeanglesect(molecule *mol, unsigned molt, unsigned moli, const char *fstr);


// Dihedral information
void readtop_dihedrals(molecule *mol, const char *fstr, unsigned molt);
void writedihedralheader(const char *fstr);
void writedihedralsect(molecule *mol, unsigned molt, unsigned moli, const char *fstr);


/*********************** MAIN ******************************/


void sep_sfg(int argc, char **argv){

  if ( (argc-3)%3 != 0 || argc < 6 )
    merror("Usage: sep_sfg <xyz file> <top file> <ntypes> ... <dens> <seed>",
	   __func__, __LINE__);

  // Number of different molecule type
  int ntypes = (argc-3)/3;

  // Getting the number of molecules of each type and the total
  // number of molecules
  unsigned nmol_ptype[ntypes];   int nmol_tot = 0;
  for ( int n=0; n<ntypes; n++ ) {
    nmol_ptype[n] = atoi( argv[n*3 + 3] ); 
    nmol_tot += nmol_ptype[n]; 
  }
  
  printpre(nmol_tot, ntypes);

  // Read the information and calculate the number of 
  // atoms in the system 
  molecule *mols = allmem(ntypes);
  unsigned natoms = 0;
  for ( int n=0; n<ntypes; n++ ){
    printinsect(n, ntypes);
    natoms += readxyz(mols, argv[n*3 + 1], n, nmol_ptype[n])*nmol_ptype[n];
    readtop_bonds(mols, argv[n*3+2], n);
    readtop_angles(mols, argv[n*3+2], n);
    readtop_dihedrals(mols, argv[n*3+2], n);
  }
  
  float dens = atof(argv[argc-2]);
  int seed = atoi(argv[argc-1]);

  // Set the lattice sites 
  lattice lat = setcubic(nmol_tot, dens);

  // Set the molecules
  unsigned *rlist = setmolpos(ntypes, nmol_tot, nmol_ptype, seed);
    
  // Write the headers in each (temporary) file 
  printf("\nWriting output files 'start.xyz' and 'start.top'..."); 
  fflush(stdout);

  writexyzheader(natoms, lat.lbox, "start.xyz");
  writebondheader("__bonds.top");
  writeangleheader("__angles.top");
  writedihedralheader("__dihedrals.top");

  // Write the config files
  for ( int n=0; n<nmol_tot; n++ ){
    writexyz(mols, rlist[n], lat, n, lat.lbox, seed, "start.xyz");
    writebondssect(mols, rlist[n], n, "__bonds.top");
    writeanglesect(mols, rlist[n], n, "__angles.top");
    writedihedralsect(mols, rlist[n], n, "__dihedrals.top");
  }
  printf("Done\n");

  // Pipe all temp. top files into start.top and remove 
  // temp. files
  pipe(mols, ntypes);

  // Free memory
  freemem(&lat, mols, ntypes, rlist);

}



/************** DEFINITIONS *****************/

void merror(const char *str, const char *fun, int line){

  fprintf(stderr, "%s at line %d, %s BAILING OUT\n", fun, line, str);

  exit(EXIT_FAILURE);

}


molecule *allmem(unsigned nmoltypes){
  
  molecule *mol = malloc(sizeof(molecule)*nmoltypes);
  if ( mol == NULL ) merror("Couldn't allocate memory",  __func__, __LINE__);
  
  return mol;
}



unsigned readxyz(molecule *mol, const char *file, 
		 unsigned moli, unsigned nmol){
  char str[256];
  float x, y, z, charge, mass;
  char type;

  FILE *fin = fopen(file, "r");
  if ( fin == NULL ) merror("Couldn't open file",  __func__, __LINE__);

  if ( fscanf(fin, "%u\n", &(mol[moli].nuau)) != 1 )
    merror("Reading failure ", __func__, __LINE__);
  
  mol[moli].rel_cm = malloc(sizeof(float)*mol[moli].nuau*3);
  mol[moli].types =  malloc(sizeof(char)*mol[moli].nuau);
  mol[moli].mass =  malloc(sizeof(float)*mol[moli].nuau);
  mol[moli].z =  malloc(sizeof(float)*mol[moli].nuau);

  if ( mol[moli].rel_cm == NULL || mol[moli].types == NULL ||
       mol[moli].z == NULL || mol[moli].mass == NULL ) 
    merror("Memory allocation error",  __func__, __LINE__);
  
  if ( fgets(str, 256, fin)== NULL ){
    merror("Couldn't read xyz-file", __func__, __LINE__);
  }
  
  float cm[3] = {0.0};
  for ( unsigned n=0; n<mol[moli].nuau; n++ ){
    
    if ( fscanf(fin, "%c%f%f%f%f%f\n", &type, &x, &y, &z, 
		 &mass, &charge) != 6 )
      merror("Reading failure ",  __func__, __LINE__);

    mol[moli].types[n] = type;
    mol[moli].rel_cm[3*n]   = x;
    mol[moli].rel_cm[3*n+1] = y;
    mol[moli].rel_cm[3*n+2] = z;

    cm[0] += x; cm[1] += y; cm[2] += z; 

    mol[moli].mass[n] = mass;
    mol[moli].z[n] = charge;
  }
  
  fclose(fin);
  
  for ( int n=0; n<3; n++ ) cm[n] /= mol[moli].nuau;

  for ( unsigned n=0; n<mol[moli].nuau; n++ ) {
    mol[moli].rel_cm[3*n]   -= cm[0];
    mol[moli].rel_cm[3*n+1] -= cm[1];
    mol[moli].rel_cm[3*n+2] -= cm[2];
  }  
  
  mol[moli].nmol_type = nmol;

  printf("Read molecule information from %s:  %d molecule(s) with %d atom(s). \n", 
	 file, mol[moli].nmol_type, mol[moli].nuau);
 
  return mol[moli].nuau;
} 


lattice setcubic(unsigned nmol_tot, float dens){
  lattice lat;

  int ngdim = ceil( pow(nmol_tot, 1.0/3.0) );
  int ngrid = ngdim*ngdim*ngdim;

  lat.z = malloc(sizeof(float)*ngrid);
  lat.y = malloc(sizeof(float)*ngrid);
  lat.x = malloc(sizeof(float)*ngrid);  

  if ( lat.z == NULL || lat.y == NULL || lat.x == NULL )
    merror("Allocation failure",  __func__, __LINE__);

  float lbox = pow( nmol_tot/dens, 1.0/3.0 );  
  float dr = lbox/ngdim;
  float hdr = 0.5*dr;

  int i = 0;
  for ( int n=0; n<ngdim; n++ ){
    float z = dr*n + hdr;
    for ( int m=0; m<ngdim; m++ ) {
      float y = dr*m + hdr;
      for ( int k=0; k<ngdim; k++ ) {
	float x = dr*k + hdr;

	lat.x[i] = x; lat.y[i] = y; lat.z[i] = z;
	i ++;
      }
    }
  }

  lat.lbox = lbox;

  return lat;
}


void writexyz(molecule *mol, unsigned moli,  lattice lat, unsigned  lati, 
	      float lbox,  int seed, const char *str){
  FILE *fptr = fopen(str, "a");
  if ( fptr == NULL ) merror("Couldn't open file",  __func__, __LINE__);
  
  const unsigned i = moli;

  seed += 10;
  srand(seed);

  for ( unsigned n=0; n<mol[i].nuau; n++ ){
    float x = lat.x[lati] + mol[i].rel_cm[3*n];  
    x=wrap(x, lbox);
    float y = lat.y[lati] + mol[i].rel_cm[3*n+1];  
    y=wrap(y, lbox);
    float z = lat.z[lati] + mol[i].rel_cm[3*n+2];  
    z=wrap(z, lbox);
     
    float vx = RAND()-0.5;
    float vy = RAND()-0.5;
    float vz = RAND()-0.5;

    fprintf(fptr, "%c %f %f %f %f %f %f %f %f\n",
	    mol[i].types[n], 
	    x, y, z,
	    vx, vy, vz,
	    mol[i].mass[n], mol[i].z[n]);
  }

  fclose(fptr);
}

unsigned *setmolpos( unsigned ntypes, unsigned nmol, unsigned *ntypes_per_mol, 
		     unsigned seed ){

  unsigned *rlist = malloc(sizeof(unsigned)*nmol);
  if ( rlist == NULL ) 
    merror("Couldn't allocate memory", __func__, __LINE__);

  unsigned i=0;
  for ( unsigned n=0; n<ntypes; n++ ){
    for ( unsigned m=0; m<ntypes_per_mol[n]; m++ ){
      rlist[i] = n;
      i++;
    } 
  }
  
  // We can mix the molecules
  if ( seed > 0 ){
    srand(seed);
    for ( unsigned n=0; n<10*nmol; n++ ){
      int a = nmol*(double)rand()/(RAND_MAX-1);
      int b = nmol*(double)rand()/(RAND_MAX-1);

      unsigned tmp = rlist[a];
      rlist[a] = rlist[b];
      rlist[b] = tmp;
    }
  }
    
  return rlist;
}

void freemem(lattice *lptr, molecule *mptr, unsigned nmols, 
	     unsigned *rlist){

  free(lptr->x); free(lptr->y); free(lptr->z); 
  
  for ( unsigned n=0; n<nmols; n++ ){
    free(mptr[n].types); 
    free(mptr[n].rel_cm);

    free(mptr[n].mass); free(mptr[n].z);
   
    free(mptr[n].bonds);
    free(mptr[n].angles);
  }
 
  free(rlist);
}


void writexyzheader(unsigned natoms, float lbox, const char *str){
  
  
  FILE *fptr = fopen(str, "w");
  if ( fptr == NULL )
    merror(" Couldn't open output file",  __func__, __LINE__);
  
  fprintf(fptr, "%u\n", natoms);
  fprintf(fptr, "%f %f %f \n", lbox, lbox, lbox);

  fclose(fptr);

}

FILE *setfilepointer(FILE *fptr, const char *section){
  char line[256];
  
  do {

    if ( fgets(line, 256, fptr) == NULL ) break;

    if ( strcmp(line, section) == 0 ){
      if ( fgets(line, 256, fptr)==NULL )
	merror("Something wrong with the reading",  __func__, __LINE__);
      return fptr;
    }

  }  while ( !feof(fptr) );

  return NULL;
}


void readtop_bonds(molecule *mol, const char *fstr, unsigned molt){
  char section[] =  "[ bonds ]\n";
  char line[256];
  fpos_t pos_file;
  unsigned a, b, btype, dummy;

  mol[molt].nbonds = 0;
  mol[molt].bonds = malloc(3*sizeof(unsigned)*mol[molt].nbonds);
  if ( mol[molt].bonds == NULL )
    merror("Memory allocation error", __func__, __LINE__);
  
  FILE *fptr = fopen(fstr, "r");
  if ( fptr == NULL ) merror("Couldn't open file",  __func__, __LINE__);

  FILE *tmpfptr = setfilepointer(fptr, section);

  if ( tmpfptr != NULL ){
    
    fptr = tmpfptr;

    do {
      fgetpos(fptr, &pos_file); 

      if ( fgets(line, 256, fptr) == NULL  ) 
	merror("Reading failure",  __func__, __LINE__);

      else if ( line[0] == '[' || line[0] == 'A' )
	break;

      else {
	fsetpos(fptr, &pos_file); 
      
	int sc = fscanf(fptr, "%u%u%u%u\n", &dummy, &a, &b, &btype);
	if ( sc != 4 ) merror("Reading failure",  __func__, __LINE__);


	(mol[molt].nbonds)++;

	mol[molt].bonds = realloc(mol[molt].bonds, 
				  3*sizeof(unsigned)*mol[molt].nbonds);
	if ( mol[molt].bonds == NULL ) 
	  merror("Memory allocation failure",  __func__, __LINE__);
	
	unsigned n = mol[molt].nbonds-1;
	mol[molt].bonds[3*n] = a; 
	mol[molt].bonds[3*n+1] = b; 
	mol[molt].bonds[3*n+2] = btype;
      
      }
    } while ( !feof(fptr) );
  }

  fclose(fptr);
   
  printf("Read bond information from %s: Found %d bond(s)\n", 
	 fstr, mol[molt].nbonds);
 
}


void writebondheader(const char *fstr){
  
  FILE *fptr = fopen(fstr, "w");
  if ( fptr == NULL ) merror("Couldn't open file", __func__, __LINE__);
  
  fprintf(fptr, "[ bonds ]\n");
  fprintf(fptr, ";WARNING generated with sfg - NO WARRANTY\n");

  fclose(fptr);
}

void writebondssect(molecule *mol, unsigned molt, unsigned moli, 
		    const char *fstr){

  static unsigned atomi = 0; //ACHTUNG!!! 

  FILE *fptr = fopen(fstr, "a");
  if ( fptr == NULL ) merror("Couldn't file", __func__, __LINE__);
 
  for ( unsigned n=0; n<mol[molt].nbonds; n++ ){
    unsigned a = mol[molt].bonds[3*n] + atomi;
    unsigned b = mol[molt].bonds[3*n+1]+ atomi;
    unsigned btype = mol[molt].bonds[3*n+2];

    fprintf(fptr, "%u %u %u %u\n", moli, a, b, btype);
  }

  fclose(fptr);

  atomi += mol[molt].nuau;
}



void readtop_angles(molecule *mol, const char *fstr, unsigned molt){
  char section[] =  "[ angles ]\n";
  char line[256];
  fpos_t pos_file;
  unsigned a, b, c, atype, dummy;

  mol[molt].nangles = 0;
  mol[molt].angles = malloc(4*sizeof(unsigned)*mol[molt].nangles);

  FILE *fptr = fopen(fstr, "r");
  if ( fptr == NULL ) merror("Couldn't open file", __func__, __LINE__);

  FILE *tmpfptr = setfilepointer(fptr, section);

  if ( tmpfptr != NULL ){
    
    fptr = tmpfptr;

    do {
      fgetpos(fptr, &pos_file); 

      if ( fgets(line, 256, fptr) == NULL  ) 
	merror("Reading failure",  __func__, __LINE__);

      else if ( line[0] == '[' )
	break;

      else {
	fsetpos(fptr, &pos_file); 
      
	int sc = fscanf(fptr, "%u%u%u%u%u\n", &dummy, &a, &b, &c, &atype);
	if ( sc != 5 ) merror("Reading failure", __func__, __LINE__);

	(mol[molt].nangles)++;

	mol[molt].angles = realloc(mol[molt].angles, 
				  4*sizeof(unsigned)*mol[molt].nangles);
	if ( mol[molt].angles == NULL ) 
	  merror("Memory allocation failure", __func__, __LINE__);
	
	unsigned n = mol[molt].nangles-1;
	mol[molt].angles[4*n] = a; 
	mol[molt].angles[4*n+1] = b;
 	mol[molt].angles[4*n+2] = c;
	mol[molt].angles[4*n+3] = atype;
      
      }
    } while ( !feof(fptr) );
  }

  printf("Read angle information from %s: Found %d angles(s)\n", 
	 fstr, mol[molt].nangles);

  fclose(fptr);
    
}




void writeangleheader(const char *fstr){
  
  FILE *fptr = fopen(fstr, "w");
  if ( fptr == NULL ) merror("Couldn't open file", __func__, __LINE__);
  
  fprintf(fptr, "[ angles ]\n");
  fprintf(fptr, ";WARNING generated with sfg - NO WARRANTY\n");

  fclose(fptr);
}


void writeanglesect(molecule *mol, unsigned molt, unsigned moli, const char *fstr){

  static unsigned atomi = 0; //ACHTUNG!!! 

  FILE *fptr = fopen(fstr, "a");
  if ( fptr == NULL ) merror("Couldn't file", __func__, __LINE__);
  
  for ( unsigned n=0; n<mol[molt].nangles; n++ ){
    unsigned a = mol[molt].angles[4*n] + atomi;
    unsigned b = mol[molt].angles[4*n+1]+ atomi;
    unsigned c = mol[molt].angles[4*n+2]+ atomi;
    unsigned atype = mol[molt].angles[4*n+3];

    fprintf(fptr, "%u %u %u %u %u\n", moli, a, b, c, atype);
  }

  fclose(fptr);

  atomi += mol[molt].nuau;
}



void readtop_dihedrals(molecule *mol, const char *fstr, unsigned molt){
  char section[] =  "[ dihedrals ]\n";
  char line[256];
  fpos_t pos_file;
  unsigned a, b, c, d, dtype, dummy;

  mol[molt].ndihedrals = 0;
  mol[molt].dihedrals = malloc(5*sizeof(unsigned)*mol[molt].ndihedrals);

  FILE *fptr = fopen(fstr, "r");
  if ( fptr == NULL ) merror("Couldn't open file",  __func__, __LINE__);

  FILE *tmpfptr = setfilepointer(fptr, section);

  if ( tmpfptr != NULL ){
    
    fptr = tmpfptr;
    
    do {
      fgetpos(fptr, &pos_file); 

      if ( fgets(line, 256, fptr) == NULL  ) 
	merror("Reading failure", __func__, __LINE__);

      else if ( line[0] == '[' )
	break;

      else {
	fsetpos(fptr, &pos_file); 
      
	int sc = fscanf(fptr, "%u%u%u%u%u%u\n", &dummy, &a, &b, &c, &d, &dtype);
	if ( sc != 6 ) merror("Reading failure", __func__, __LINE__);


	(mol[molt].ndihedrals)++;

	mol[molt].dihedrals = realloc(mol[molt].dihedrals, 
				  5*sizeof(unsigned)*mol[molt].ndihedrals);
	if ( mol[molt].dihedrals == NULL ) 
	  merror("Memory allocation failure", __func__, __LINE__);
	
	unsigned n = mol[molt].ndihedrals-1;
	mol[molt].dihedrals[5*n] = a; 
	mol[molt].dihedrals[5*n+1] = b; 
	mol[molt].dihedrals[5*n+2] = c; 
	mol[molt].dihedrals[5*n+3] = d; 
	mol[molt].dihedrals[5*n+4] = dtype;
      
      }
    } while ( !feof(fptr) );
  }

  fclose(fptr);
   
  printf("Read dihedrals information from %s: Found %d dihedral(s)\n", 
	 fstr, mol[molt].ndihedrals);
 
}


void writedihedralheader(const char *fstr){
  
  FILE *fptr = fopen(fstr, "w");
  if ( fptr == NULL ) merror("Couldn't open file", __func__, __LINE__);
  
  fprintf(fptr, "[ dihedrals ]\n");
  fprintf(fptr, ";WARNING generated with sfg - NO WARRANTY\n");

  fclose(fptr);
}

void writedihedralsect(molecule *mol, unsigned molt, unsigned moli, const char *fstr){

  static unsigned atomi = 0; //ACHTUNG!!! 

  FILE *fptr = fopen(fstr, "a");
  if ( fptr == NULL ) merror("Couldn't file", __func__, __LINE__);
 
  for ( unsigned n=0; n<mol[molt].ndihedrals; n++ ){
    unsigned a = mol[molt].dihedrals[5*n] + atomi;
    unsigned b = mol[molt].dihedrals[5*n+1] + atomi;
    unsigned c = mol[molt].dihedrals[5*n+2] + atomi;
    unsigned d = mol[molt].dihedrals[5*n+3] + atomi;
    unsigned dtype = mol[molt].dihedrals[5*n+4];

    fprintf(fptr, "%u %u %u %u %u %u\n", moli, a, b, c, d, dtype);
  }

  fclose(fptr);

  atomi += mol[molt].nuau;
}

void pipe( molecule *mol, int ntypes ){

  int bonds = 0;  int ang = 0;  int dihed = 0;
  for ( int n=0; n<ntypes; n++ ){
    bonds += mol[n].nbonds;
    ang += mol[n].nangles;
    dihed += mol[n].ndihedrals;
  }

  if ( system("rm -f start.top") == -1 )
    merror("Linux system call 'rm' failed", __func__, __LINE__);

  int retsys = 0;
  if ( bonds > 0 ){
    retsys = system("cat __bonds.top >> start.top");
    if ( retsys == -1 )
      merror("Linux system call 'rm' failed", __func__, __LINE__);

    retsys = system("echo ""  >> start.top");
    if ( retsys == -1 )
      merror("Linux system call 'rm' failed", __func__, __LINE__);
    
  }
  if ( ang > 0 ) {
    retsys = system("cat __angles.top >> start.top");
    if ( retsys == -1 )
      merror("Linux system call 'rm' failed", __func__, __LINE__);

    retsys = system("echo ""  >> start.top");
    if ( retsys == -1 )
      merror("Linux system call 'rm' failed", __func__, __LINE__);

  }
  if ( dihed > 0 ){
    retsys=system("cat __dihedrals.top >> start.top");
    if ( retsys == -1 )
      merror("Linux system call 'rm' failed", __func__, __LINE__);

    retsys = system("echo ""  >> start.top");
    if ( retsys == -1 )
      merror("Linux system call 'rm' failed", __func__, __LINE__);
  }

}


float velocity(void){

  return (float)rand()/(RAND_MAX - 1) - 0.5;
    
}


int findnumtypes(molecule *mol, int ntypes_mol){

  int type = -1;
  for ( int n=0; n<ntypes_mol; n++ ){
    for ( unsigned m=0; m<mol[n].nuau; m++ )
      if ( (signed)mol[n].types[m] > type ) type = mol[n].types[m];
  }

  return type + 1;

}

void printpre(int nmol_tot, int ntypes ){

  printf("This is sep_sfg (System File Generator) version %s\n", VERSION);
  printf("This program is a part of the sep library - there is absolutely ");
  printf("NO WARRANTY\n\n");

  
  printf("From arguments: Number of molecules: %d [%d type(s)]\n", 
	 nmol_tot, ntypes);

}

void  printinsect(int n, int ntypes){
  
  printf("\n------------[info] Type %d of %d------------\n", n+1, ntypes);

}

float wrap(float x, float lbox){

  if ( x>lbox ) x -= lbox;
  else if ( x<0.0 ) x += lbox;

  return x;
}
