
#ifndef __TASK_H__
#define __TASK_H__

#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include "sep.h"

typedef struct {

  int type;
  int block;

  bool pairflag; 
  double max_cutoff;
  
  // LJ 
  char ptypes[3];
  double cutoff;
  double sigma;
  double epsilon;

  // Bond 
  int btype;
  double lbond;
  double kbond;

  // Angle
  int atype;
  double angle;
  double kangle;

  // Diheadral
  int dtype;
  double dparam[6];

  // Coulomb
  double sf_cutoff;
  
} taskmanager;


int tasktoint(char *taskstr);
void printtask(taskmanager *ptask, int tasknr);
void settask(taskmanager *ptask, int block, char *taskopt, ... );
void dotask(double **f, seppart *atoms, taskmanager *ptask,
	    int tid, sepsys *sys);
void dotask2(seppart *pptr, taskmanager *ptask, int ntasks, sepsys *sys,
	      const unsigned exclopt);
void dotask3(seppart *pptr, taskmanager *ptask, int ntasks, sepsys *sys,
	     const unsigned exclopt);
void dotask4(seppart *pptr, taskmanager *ptask, int ntasks, sepsys *sys,
	     const unsigned exclopt);
// Octave interface
// molsim('task', 'lj', 'AA', 2.5, 1.0, 1.0, 1, 1);
// -> block , total # tasks 
// molsim('task', 'bond', 0, 1.0, 500, 2, 3);
// molsim('task', 'angle', 0, pi, 100, 2, 3);

#endif
