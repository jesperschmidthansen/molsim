
/* 
* sepdef.h - This file is a part of the sep-library 
*
* Copyright (C) 2008 Jesper Schmidt Hansen 
* 
* License: GPL - see COPYING for copying conditions.
* There is ABSOLUTELY NO WARRANTY, not even for MERCHANTIBILITY or
* FITNESS FOR A PARTICULAR PURPOSE.
*
* Contact: schmidt@zigzak.net
*/

#ifndef __SEPDEF_H__
#define __SEPDEF_H__

#define SEP_FALSE 0
#define SEP_TRUE 1

#define SEP_BOND 10

#define SEP_NEIGHB 3000
#define SEP_NUM_NEIGHB 3000

#define SEP_NO_NEIGHB 0

#define SEP_NEIGHB_ALL 1
#define SEP_ALL 1

#define SEP_NEIGHB_EXCL_BONDED 2
#define SEP_EXCL_BONDED 2
        
#define SEP_NEIGHB_EXCL_SAME_MOL 3
#define SEP_EXCL_SAME_MOL 3
        
#define SEP_MAX_NUM_MOL 10000

#define SEP_PI 3.14159265358979

#define SEP_WCACF 1.12246204830937 

#define SEP_LJCF2 0.016316891136 

#define SEP_R250_LEN     250
#define SEP_R521_LEN     521

#define SEP_MAXNID 100000 
#define SEP_NRETVALS 10

#define SEP_BRUTE 0
#define SEP_NEIGHBLIST 1
#define SEP_LLIST_NEIGHBLIST 2

#define SEP_LEAPFROG 0
#define SEP_NOSEHOOVER 1
#define SEP_ANDERSEN 2
#define SEP_SHAKE 3

#define SEP_SHAKE_MAXIT 1000

#define SEP_SUCCESS 1
#define SEP_FAILURE 1

#ifdef COMPLEX

#include <complex.h>

#ifdef imaginary
#define SEP_I _Imaginary_I
#else
#define SEP_I _Complex_I
#endif

#endif



#endif
