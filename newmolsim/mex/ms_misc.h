#ifndef __MISC__
#define __MISC__

#define _Wrap( x, y )              \
{                                  \
if ( x > 0.5*y ) x -= y;           \
else if  ( x < -0.5*y ) x += y;    \
}



#endif
