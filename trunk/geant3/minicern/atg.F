*
* $Id: atg.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: atg.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:44  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:33  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/04/01 15:01:49  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION ATG(X1,X2)

      PARAMETER (PI = 3.14159 26535 89793D0)

      ATG=ATAN2(X1,X2)+(PI-SIGN(PI,X1))
      RETURN
      END
