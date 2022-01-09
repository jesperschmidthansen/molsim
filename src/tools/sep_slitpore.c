/***************************************************************
 * sep_slitpore.c
 *
 * Sets up a slitpore xyz file. Atomic systems only. The wall is 
 * simple cubic lattice; three wall layers. 
 *
 * usage: sep_slitpore <nxyz> <Wall density> <Fluid density>
 *
 * Author: Jesper Hansen
 ***************************************************************/

#include "sep.h"

int main(int argc, char **argv){
  sepatom *atoms;
  sepsys sys;

  
  if ( argc != 5 )
    sep_error("Usage: sep_slitpore <nxyz> <Wall density> <Fluid density> <opt (1 for two walls otherwise one wall)>");

  srand(21212);
 
  // Setting parameter values
  int nxyz = atoi(argv[1]);
  double densW = atof(argv[2]); 
  double densF = atof(argv[3]);
  int opt = atof(argv[4]);
 
  if ( densF > densW )
    sep_error("Fluid density must be smaller than wall density - sorry");
  
  int natoms = nxyz*nxyz*nxyz;
  double lbox = pow(natoms/densW, 1.0/3.0);
  
  // Allocating memory 
  atoms = sep_init(natoms, SEP_NEIGHB);

  // Setting up the system
  sys = sep_sys_setup(lbox, lbox, lbox, 2.5, 0.0, natoms, SEP_BRUTE);
  
  // Initializing the positions and momenta
  sep_set_lattice(atoms, sys);

  int nW = nxyz*nxyz*3;
   
  for ( int n=0; n<nW; n++ )
    atoms[n].type = 'W';

  if ( opt == 1 )  {
    for ( int n=sys.npart-1; n>sys.npart-nW-1; n-- )
      atoms[n].type = 'w';
  }

  int nF;
  if ( opt == 1 )
    nF = (natoms - 2*nW)*densF/densW;
  else
    nF = (natoms - nW)*densF/densW;
  
  int cF = 0;
  while ( cF < nF ){
    int index = ceil(sep_rand()*natoms)-1;
    if ( atoms[index].type == 'A' ){
      atoms[index].type = 'F';
      cF ++;
    }
  }

  sep_save_xyz(atoms, "wWF", "slitpore.xyz", "w", sys);

  // Freeing memory 
  sep_close(atoms, natoms);
  sep_free_sys(&sys);

  return 0;
} 
