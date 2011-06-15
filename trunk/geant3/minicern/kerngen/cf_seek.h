/*
* $Id: cf_seek.h 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cf_seek.h,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.2  1997/02/04 17:35:36  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:26  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:17  mclareni
* Kernlib
*
*
*
* cf#seek.inc
*/
#if defined(CERNLIB_QMAPO)
#include <sys/file.h>        /*  Apollo                     */
#elif defined(CERNLIB_QMAMX)
#include <sys/types.h>       /*  AMX                        */
#include <sys/file.h>
#elif defined(CERNLIB_QMOS9)
#include <stdio.h>           /*  Microware OS-9             */
#elif defined(CERNLIB_QMVAX)
#include <file.h>            /*  VAX/VMS                    */
#elif defined(CERNLIB_QMIRTD)
#include <unistd.h>          /*  IRTD                */
#elif defined(CERNLIB_QMDOS) || defined(CERNLIB_WINNT)
 #ifdef __GNUC__
  #include <sys/file.h>
  #include <unistd.h>
 #else
  #ifdef __STDC__
   #undef __STDC__
  #endif
  #include <stdio.h>
  #ifdef WIN32
   #include <io.h>
  #endif
 #endif
#else
#include <sys/types.h>       /*  default Posix              */
#include <unistd.h>
#endif
