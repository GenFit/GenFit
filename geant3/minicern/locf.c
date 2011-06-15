/*
 * $Id: locf.c 220 2007-11-19 16:08:06Z rdm $
 *
 * $Log: locf.c,v $
 * Revision 1.1.1.1  2002/07/24 15:56:28  rdm
 * initial import into CVS
 *
 * Revision 1.1.1.1  2002/06/16 15:18:46  hristov
 * Separate distribution  of Geant3
 *
 * Revision 1.1.1.1  1999/05/18 15:55:28  fca
 * AliRoot sources
 *
 * Revision 1.3  1997/09/02 14:26:38  mclareni
 * WINNT correction
 *
 * Revision 1.2  1997/02/04 17:34:35  mclareni
 * Merge Winnt and 97a versions
 *
 * Revision 1.1.1.1.2.1  1997/01/21 11:29:36  mclareni
 * All mods for Winnt 96a on winnt branch
 *
 * Revision 1.1.1.1  1996/02/15 17:49:24  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
#include "kerngen/fortranc.h"
#if defined(CERNLIB_LXIA64)
#include "stdio.h"
#endif

#if defined(CERNLIB_MSSTDCALL) && defined(CERNLIB_LOCF_CHARACTER)
# define Dummy2LocPar  ,_dummy
# define DummyDef     int _dummy;
#else
# define Dummy2LocPar  
# define DummyDef
#endif

#if defined(CERNLIB_QMIRTD)
#include "irtdgs/locf.c"
#elif defined(CERNLIB_QMVAOS)
#include "vaogs/locf.c"
#else
/*>    ROUTINE LOCF
  CERN PROGLIB# N100    LOCF            .VERSION KERNFOR  4.36  930602
*/
#define NADUPW 4   /* Number of ADdress Units Per Word */
#define LADUPW 2   /* Logarithm base 2 of ADdress Units Per Word */
#if defined(CERNLIB_QX_SC)
unsigned int type_of_call locf_(iadr Dummy2LocPar)
#elif defined(CERNLIB_QXNO_SC)
unsigned int type_of_call locf(iadr Dummy2LocPar)
#elif defined(CERNLIB_QXCAPT)
unsigned int type_of_call LOCF(iadr Dummy2LocPar)
#endif
   char *iadr;
#ifdef DummDef
   DummyDef
#endif
{
#if defined(CERNLIB_LXIA64)
  const unsigned long long int mask=0x00000000ffffffff;
  static unsigned long long int base=1;
  unsigned long long int jadr=(unsigned long long int) iadr;
  unsigned long long int jadrl = ((mask & jadr) >> LADUPW);

  if (base == 1) {
    base = (~mask & jadr);
  } else if(base != (~mask & jadr)) {
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("locf_() Warning: changing base from %lx to %lx!!!\n",
    	   base, (~mask & jadr));
    printf("This may result in program crash or incorrect results\n");
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  }
  return ((unsigned) jadrl);
#else
  return( ((unsigned) iadr) >> LADUPW );
#endif
}
#undef Dummy2LocPar
#undef DummyDef
#undef CERNLIB_LOCF_CHARACTER
/*> END <----------------------------------------------------------*/
#endif

