/*
 * $Id: ishftr.c 220 2007-11-19 16:08:06Z rdm $
 *
 * $Log: ishftr.c,v $
 * Revision 1.1.1.1  2002/06/16 15:18:47  hristov
 * Separate distribution  of Geant3
 *
 * Revision 1.1.1.1  1999/05/18 15:55:33  fca
 * AliRoot sources
 *
 * Revision 1.1.1.1  1996/02/15 17:50:07  mclareni
 * Kernlib
 *
 */
#include "kerngen/pilot.h"
/*>    ROUTINE ISHFT
  CERN PROGLIB#         ISHFTR          .VERSION KERNLNX  1.02  940511

  Logical right shift by *len (+ve) places
*/
unsigned int ishftr_(arg,len)
unsigned int *arg;
int *len;
{
   return(*arg >> *len);
}
/*> END <----------------------------------------------------------*/
