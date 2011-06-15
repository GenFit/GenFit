/*
* $Id: qf_f2c.h 220 2007-11-19 16:08:06Z rdm $
*
* $Log: qf_f2c.h,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:20  mclareni
* Kernlib
*
*
*  This is a piece of f2c.h file by AT&T
*
* qf_f2c.h
*/
typedef long ftnlen;
typedef long int integer;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;

typedef struct { real r, i; } complex;
#ifdef WIN32               /*   V.E.Fine 03.07.93 */
#define _COMPLEX_DEFINED
#endif

typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
/* typedef long long longint; */ /* system-dependent */

#define TRUE_ (1)
#define FALSE_ (0)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
