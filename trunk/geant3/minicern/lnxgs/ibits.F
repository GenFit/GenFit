*
* $Id: ibits.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ibits.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:30  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      FUNCTION IBITS (IZW,IZP,NZB)
*
      IBITS=ISHFTR(LSHIFT(IZW,32-IZP-NZB),32-NZB)
*
      RETURN
      END
