

/* 
* sepret.h - This file is a part of the sep-library 
*
* Copyright (C) 2011 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/

#ifndef __SEPRET_H__
#define __SEPRET_H__


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>
#include <float.h>

#include "sepstrct.h"
#include "sepmol.h"

/** 
 * Resets the members in return structure. Called prior to the force functions
 * and integrator
 * @param retval Pointer to the seplib return structure
 */
void sep_reset_retval(sepret *retval);

/**
 * Retrieves the normal atomic pressure
 * @param retval Pointer to the seplib return structure
 * @param sys Pointer to the seplib system structure
 * @return The atomistic pressure
 */
double sep_get_pressure(sepret *retval, sepsys *sys);

/**
 * Retrieves the kinetic temperature
 * @param retval Pointer to the seplib return structure
 * @param sys Pointer to the seplib system structure
 * @return The kinetic temperature
 */
double sep_get_temperature(sepret *retval, sepsys *sys);


/**
 * Calculates the atomic pressure tensor. Called after force functions and 
 * integrators
 * @param retval Pointer to the seplib return structure
 * @param sys Pointer to the seplib system structure
 */
void sep_pressure_tensor(sepret *retval, sepsys *sys);

/**
 * Calculates the molecular pressure tensor. Called after force functions 
 * and integrators 
 */
void sep_mol_pressure_tensor(sepatom *atoms, sepmol *mols,
			     sepret *ret, sepsys *sys);


#endif
