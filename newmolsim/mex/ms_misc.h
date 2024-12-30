#ifndef __MISC__
#define __MISC__

#define MAX_NNEIGHB 500
#define MAX_EXCLUSIONS 5



#define _Wrap( x, y )	               \
{    	                               \
	if ( x > 0.5*y ) x -= y;           \
	else if  ( x < -0.5*y ) x += y;    \
}

#define _Periodic0( x, y )             \
{                                      \
	if ( x > y ) x -= y;               \
 	else if  ( x < 0.0f ) x += y;      \
}


#define _Periodic( cross, x, y )                  \
{                                                 \
 	cross = 0;                                    \
	if ( x > y ) { x -= y; cross = 1; }           \
 	else if  ( x < 0.0f ) { x += y; cross = -1;}  \
}

#define _Dot3( c, a, b )  \
	c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] 

#endif
