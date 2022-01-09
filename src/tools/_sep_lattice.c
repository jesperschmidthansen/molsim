

/***************************************************************
 * sep_lattice.c
 *
 * Build a cubic or orthombic lattice and write to file start.xyz
 *
 * Copyright 2021 Jesper Hansen
 ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define DEFAULT_NXYZ 10
#define DEFAULT_LENGTH 10.0

void printhelp(void){

  printf("\nA simple orthombic lattice generator  - Copyright Jesper Hansen\n\n");
  printf("usage: sep_lattice [-nlswWbrft] \n\n");
  printf("Options\n");
  printf("-n: Number of particles in each direction e.g. -n=12,14,9 (default 10) ");
  printf("If one number is specified this is used in all directions\n");
  printf("-l: Box length in each direction e.g. -l=10,21.1,15.511 (default 10.0) ");
  printf("If one number is specified this is used in all directions\n");
  printf("-s: Seed for velocities e.g. -s=42 (default time(NULL))\n");
  printf("-w: Option for single wall e.g. -w=0.6 (number is fluid density)\n");
  printf("-W: Option for double wall e.g. -W=0.6 (number is fluid density)\n");
  printf("-b: Option for body centered lattice (default false)\n");
  printf("-r: Specifying the density (leaves -l option void and uses first value for -n option)\n");
  printf("-f: Output file (default start.xyz)\n");
  printf("-t: Set type specifications. ntypes,labels,numbers,masses, e.g. -t=2,A,B,800,200,1.0,3.0\n");
  printf("\n");
  
  printf("Examples\n");
  printf("sep_lattice -n=10,10,10 -l=10.0,10.0,10.0 (default and same as sep_lattice)\n");	 
  printf("sep_lattice -n=25,25,25 -l=25,25,30 -W=0.75 -b -s=42 -f=mystartfile.xyz\n");	 
    
}


double getdensity(const char *str){
  double density; 

  sscanf(str, "-r=%lf", &density);

  return density;
}

unsigned gettypes(char *types, unsigned *ntypes, double *mass, const char *str){
  unsigned int numtype;
  
  sscanf(str, "-t=%u", &numtype);

  if ( numtype == 2 ){
   int nread = sscanf(str, "-t=%u,%c,%c,%u,%u,%lf,%lf", &numtype,
		      &types[0], &types[1], &ntypes[0], &ntypes[1],
		      &mass[0], &mass[1]);

   if ( nread != 7 ){
     fprintf(stderr, "Error reading input from -t option\n");
     exit(EXIT_FAILURE);
   }
   
   types[2]='\0';
  }
  else if ( numtype == 3 ){
    int nread = sscanf(str, "-t=%u,%c,%c,%c,%u,%u,%u,%lf,%lf,%lf", &numtype,
		       &types[0], &types[1], &types[2],
		       &ntypes[0], &ntypes[1], &ntypes[3],
		       &mass[0], &mass[1], &mass[2]);

    if ( nread != 10 ){
      fprintf(stderr, "Error reading input from -t option\n");
      exit(EXIT_FAILURE);
    }
   
    types[3]='\0';
  }
  else {
    fprintf(stderr, "Sorry only up to three types are supported\n");
    exit(EXIT_FAILURE);
  }

  return numtype;
}


void getlengths(double *lbox, const char *str){

  bool flag = false;

  for ( size_t n=0; n<strlen(str); n++ ){
    if ( str[n] == ',' ) {
      flag = true;
      break;
    }
  }
      
  if ( !flag ){
    sscanf(str, "-l=%lf", &lbox[0]);
    lbox[1]=lbox[2]=lbox[0];
  }
  else {
    sscanf(str, "-l=%lf,%lf,%lf", &lbox[0], &lbox[1], &lbox[2]);
  }
  
}

void getnxyz(int *nparts, const char *str){
  bool flag = false;

  for ( size_t n=0; n<strlen(str); n++ ){
    if ( str[n] == ',' ) {
      flag = true;
      break;
    }
  }

  if ( !flag ){
    sscanf(str, "-n=%d", &nparts[0]);
    nparts[1]=nparts[2]=nparts[0];
  }
  else {
     sscanf(str, "-n=%d,%d,%d", &nparts[0], &nparts[1], &nparts[2]);
  }
  
}


void sep_lattice( int argc, char **argv ){
  int nxyz[3], mseed = time(NULL);
  unsigned numtype=0, ntypes[3];
  double lbox[3], fdens = 100.0, mass[3];
  bool wflag = false;
  bool Wflag = false;
  bool bcflag = false;
  bool tflag = false;
  char filestr[256] = "start.xyz";
  char types[4];

  // Default values
  nxyz[0] = nxyz[1] = nxyz[2] = DEFAULT_NXYZ;
  lbox[0] = lbox[1] = lbox[2] = DEFAULT_LENGTH;

  for ( int n=1; n<argc; n++ )
    printf("%s\n", argv[n]);
      
  // Get input values
  if ( argc > 1 ){

    for ( int n=1; n<argc; n++ ){

      char input = argv[n][1];
   
      if ( input == 'n' )
	getnxyz(nxyz, argv[n]);
      else if ( input == 'l' )
	getlengths(lbox, argv[n]);
      else if ( input == 'r' ){
	double density = getdensity(argv[n]);
	nxyz[1]=nxyz[2]=nxyz[0];
	lbox[0]=lbox[1]=lbox[2]=pow(nxyz[0]*nxyz[1]*nxyz[2]/density, 1./3.);
      }
      else if ( input == 'f' )
	sscanf(argv[n], "-f=%s", filestr);
      else if ( input == 's' )
	sscanf(argv[n], "-s=%d", &mseed);
      else if ( input == 'w' ){
	wflag = true;
	sscanf(argv[n], "-w=%lf", &fdens);
      }
      else if ( input == 'W' ){
	Wflag = true;
	sscanf(argv[n], "-W=%lf", &fdens);
      }
      else if ( input == 'b' )
	bcflag = true;
      else if ( input == 't' ){
	numtype = gettypes(types, ntypes, mass, argv[n]);
	tflag = true;
      }
      else if  ( input == 'h' ){
	printhelp();
	exit(EXIT_SUCCESS);
      }
      else {
	fprintf(stderr, "Error - use -h option for help");
	exit(EXIT_FAILURE);
      }
      
    }
  }

  unsigned long npart = nxyz[0]*nxyz[1]*nxyz[2];
   
  // Just checking for body centered (dimensions must be even)
  if ( bcflag ){
    if ( nxyz[1]%2!=0 || nxyz[2]%2!=0 ){
      fprintf(stderr, "Body centered lattice demands even particles dimensions\n");
      exit(EXIT_FAILURE);
    }
  }

  // Setting particle velocities
  double v[npart][3];
  double sumv[3]={0.0};

  srand(mseed);
  for ( long unsigned n=0; n<npart; n++ ){
    for ( int k=0; k<3; k++ ){
      v[n][k] = (rand()/(RAND_MAX+1.0)-0.5)*2.0; 
      sumv[k] += v[n][k];
    }
  }

  for ( int k=0; k<3; k++ ) sumv[k] = sumv[k]/npart;

  for ( long unsigned n=0; n<npart; n++ )
    for ( int k=0; k<3; k++ ) v[n][k] = v[n][k] - sumv[k];
  
  // Setting positions
  double gab[3], x, y, z;
  double r[npart][3];
  double pmass[npart];
  
  for ( int n=0; n<3; n++) gab[n] = lbox[n]/nxyz[n];

  char *ptype = malloc(sizeof(char)*npart); 
  if ( ptype==NULL ) {
    fprintf(stderr, "Error allocating memory");
    exit(EXIT_FAILURE);
  }
  
  long unsigned pindex = 0;
  for ( int nz=0; nz<nxyz[2]; nz++ ){
    z = nz*gab[2];
    for ( int ny=0; ny<nxyz[1]; ny++ ){
      y = ny*gab[1];
      for ( int nx=0; nx<nxyz[0]; nx++ ){
	x = nx*gab[0];
	if ( bcflag ) {
	  if ( (nz%2 != 0 && ny%2==0) || (nz%2==0 && ny%2!=0 ) )
	    x = x + gab[0]*0.5;
	}
	r[pindex][0]=x; r[pindex][1]=y; r[pindex][2]=z;
	ptype[pindex]='A';
	pmass[pindex]=1.0;
	
	pindex++;
      }
    }
  }

  // For confinement the fluid A labels must be relabled  F, W and w
  unsigned long int nF = 0; unsigned nW = 0;
  if ( wflag || Wflag ){
  
    double density = nxyz[0]*nxyz[1]*nxyz[2]/(lbox[0]*lbox[1]*lbox[2]);
    if ( density <= fdens ){
      fprintf(stderr, "Error: fluid density must be smaller than wall densities");
      exit(EXIT_FAILURE);
    }
    
    nW = nxyz[0]*nxyz[1]*3;
    for ( unsigned int n=0; n<nW; n++ ) ptype[n] = 'w';
  
    if ( Wflag ){
      nF = (npart - 2*nW)*fdens/density;
      for ( long unsigned n=npart-1; n>=npart-nW; n-- ) ptype[n] = 'W';
    }
    else if ( wflag ){
      nF = (npart - nW)*fdens/density;
    }

    unsigned cF = 0;
    while ( cF < nF ){
      int index = (double)rand()*npart/RAND_MAX;

      if ( ptype[index] == 'A' ){
	ptype[index] = 'F';
	cF ++;
      }
    }
    
  }

  // Types and masses
  if ( tflag ) {

    unsigned sumtypes = 0;
    for ( unsigned n=0; n<numtype; n++ ) sumtypes += ntypes[n];
    if ( sumtypes != npart ){
      fprintf(stderr, "Particles do not add up. Check your -t option\n");
      exit(EXIT_FAILURE);
    }
      
    for ( unsigned n=0; n<numtype; n++ ){
      for ( unsigned m=0; m<ntypes[n]; m++ ){
	ptype[n*m + m] = types[n];
	pmass[n*m + m] = mass[n];
      }
    }

    // Shuffle/Mix
    for ( unsigned n=0; n<npart; n++ ){
      int ia = (int)( (rand()/(RAND_MAX+1.0))*npart );
      int ib = (int)( (rand()/(RAND_MAX+1.0))*npart );

      double tmp_type = ptype[ia];
      double tmp_mass = pmass[ia];

      ptype[ia] = ptype[ib]; pmass[ia]=pmass[ib];
      ptype[ib] = tmp_type; pmass[ib]=tmp_mass;
    }
      
  }
  
  printf("Writing file %s... ", filestr);
  fflush(stdout);

  // Print header
  FILE *fout = fopen(filestr, "w");
  if ( fout == NULL ) {
    fprintf(stderr, "Could not open file for writing...Aaargh..."); 
    exit(EXIT_FAILURE); 
  }


  if ( wflag )
    fprintf(fout, "%lu\n", nW + nF);
  else if ( Wflag )
    fprintf(fout, "%lu\n", 2*nW + nF);
  else
    fprintf(fout, "%lu\n", npart);

  fprintf(fout, "%lf %lf %lf\n", lbox[0], lbox[1], lbox[2]);
  
  // Print conf.
  for ( unsigned long int n=0; n<npart; n++ ){
    
    if ( wflag || Wflag ){
      fprintf(fout, "%c %f %f %f %f %f %f %f %f\n", ptype[n],
	      r[n][0],r[n][1],r[n][2], v[n][0], v[n][1], v[n][2],
	      pmass[n], 0.0);
    }
    else if ( !wflag && !Wflag )
      fprintf(fout, "%c %f %f %f %f %f %f %f %f\n", ptype[n],
	      r[n][0],r[n][1],r[n][2], v[n][0], v[n][1], v[n][2],
	      pmass[n], 0.0);
  }
    
  fclose(fout);
  printf("done \n");

  free(ptype);

 
}   

