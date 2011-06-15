*
* $Id: history.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: history.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.2  2002/05/13 12:40:58  hristov
* Dummy subroutines to avoid files with no code in
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
*
* *** Version 1.21 (23-04-93)
*
*   8-04-92 : define Patchy title
*   7-04-92 : mods in GPTSET
*           : mods in START,STAMAT for file definition on ibm
*  16-04-92 : Cosmetic mods in START
*  27-04-92 : Aname bug on IBM corrected, mod of .mat output
*  10-11-92 : Added transfer of colour attribut to SET output
*  04-03-93 : Added routine GWRTRE to write the tree information into
*             the material file
*  04-03-93 : Added an option to choose the number of instances
*  14-04-93 : Mods in the output of tree information
*  14-04-93 : Mods in the "number of instances" option so that it affects
*             only to divisions
*  23-04-93 : Names changes to GEANT name convention scheme.  Some
*             cleaning done to GDEXCA (Earlier SETDRS)
      SUBROUTINE HISTORY_DUMMY
      END
