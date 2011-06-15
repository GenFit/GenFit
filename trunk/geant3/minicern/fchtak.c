/*
 * $Id: fchtak.c 220 2007-11-19 16:08:06Z rdm $
 *
 * $Log: fchtak.c,v $
 * Revision 1.1.1.1  2002/06/16 15:18:46  hristov
 * Separate distribution  of Geant3
 *
 * Revision 1.1.1.1  1999/05/18 15:55:29  fca
 * AliRoot sources
 *
 * Revision 1.1.1.1  1996/02/15 17:49:40  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE FCHTAK
  CERN PROGLIB#         FCHTAK          .VERSION KERNFOR  4.31  911111
  ORIG. 22/02/91, JZ

      copy a Fortran character string
      to allocated memory zero-terminated,
      return the memory pointer
*/
#include <stdio.h>
#include "kerngen/fortchar.h"
char *fchtak(ftext,lgtext)
#if defined(CERNLIB_QMCRY)
      _fcd  ftext;
#endif
#if !defined(CERNLIB_QMCRY)
      char *ftext;
#endif
      int  lgtext;
{
      char *malloc();
      char *ptalc, *ptuse;
      char *utext;
      int  nalc;
      int  ntx, jcol;

      nalc  = lgtext + 8;
      ptalc = malloc (nalc);
      if (ptalc == NULL)     goto exit;
#if defined(CERNLIB_QMCRY)
      utext = _fcdtocp(ftext);
#endif
#if !defined(CERNLIB_QMCRY)
      utext = ftext;
#endif

      ptuse = ptalc;
      ntx   = lgtext;
      for (jcol = 0; jcol < ntx; jcol++)  *ptuse++ = *utext++;

      *ptuse = '\0';
exit: return  ptalc;
}
/*> END <----------------------------------------------------------*/
