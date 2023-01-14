
/* 
 * sepret.c - This file is a part of the sep-library 
*
* Copyright (C) 2011 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/


#include "sepret.h"

// Note: memset does not give speed-up (compiler optimization) and
// the explicit/manual reset is preferred
void sep_reset_retval(sepret *retval){
  
  retval->epot = 0;
  retval->ecoul = 0;
  retval->ekin = 0;
  retval->sumv2 = 0;
  
  for ( int k=0; k<3; k++ ){
    for ( int kk=0; kk<3; kk++ ){
      retval->pot_P[k][kk] = 0.0;
      retval->kin_P[k][kk] = 0.0;
      retval->P[k][kk] = 0.0;
      
      retval->pot_P_conservative[k][kk]=.0;
      retval->pot_P_random[k][kk]=0.0;
      retval->pot_P_dissipative[k][kk]=0.0;
      retval->pot_P_bond[k][kk]=0.0;
      
      retval->pot_P_mol[k][kk] = 0.0;
      retval->kin_P_mol[k][kk] = 0.0;
      retval->P_mol[k][kk] = 0.0;
      
      retval->pot_T_mol[k][kk] = 0.0;
      retval->kin_T_mol[k][kk] = 0.0;
      retval->T_mol[k][kk] = 0.0;
    }
  }
  
}


double sep_get_pressure(sepret *retval, sepsys *sys){
  
  sep_pressure_tensor(retval, sys);
  
  return retval->p; 
}


double sep_get_temperature(sepret *retval, sepsys *sys){
  
  return 2.0/(3.0*sys->ndof)*retval->ekin; 

}


void sep_pressure_tensor(sepret *retval, sepsys *sys){
  
  const double ivol = 1.0/sys->volume;

  for ( int k=0; k<3; k++ ){
    for ( int kk=0; kk<3; kk++ ){
      retval->P[k][kk] = (retval->kin_P[k][kk] + retval->pot_P[k][kk])*ivol; 
    }
  }
  
  retval->p=0.0;
  for ( int k=0; k<3; k++ ){
    retval->p += retval->P[k][k];
  }

  retval->p /= 3.0;

}


void sep_mol_pressure_tensor(sepatom *atoms, sepmol *mols, sepret *retval, 
			     sepsys *sys){
  
  const double ivol = 1.0/sys->volume;
  sep_eval_mol_pressure_tensor(atoms, mols, retval, sys);

  for ( int k=0; k<3; k++ )
    for ( int kk=0; kk<3; kk++ )
      retval->P_mol[k][kk] = (retval->kin_P_mol[k][kk] + 
			      retval->pot_P_mol[k][kk])*ivol; 

  retval->p_mol=0.0;
  for ( int k=0; k<3; k++ )
    retval->p_mol += retval->P_mol[k][k];

  retval->p_mol /= 3.0;

}

