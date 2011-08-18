/*
* $Id: cf_clos.h 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cf_clos.h,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.2  1997/02/04 17:35:35  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:25  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:20  mclareni
* Kernlib
*
*
*
* cf#clos.inc
*/
#if defined(CERNLIB_QMAPO)||defined(CERNLIB_QMOS9)
#elif defined(CERNLIB_QMVAX)
#include <file.h>            /*  VAX/VMS                    */
#elif defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
  #ifdef WIN32
#   ifdef __STDC__
#    undef __STDC__
#   endif
    #include <io.h>
  #endif
#else
#include <unistd.h>          /*  default Posix              */
#endif
