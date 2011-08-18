/*
 * $Id: qnexte.c 220 2007-11-19 16:08:06Z rdm $
 *
 * $Log: qnexte.c,v $
 * Revision 1.1.1.1  1999/05/18 15:55:28  fca
 * AliRoot sources
 *
 * Revision 1.1.1.1  1996/02/15 17:49:35  mclareni
 * Kernlib
 *
 */
/*
* CERN PROGLIB# Z041    QNEXTE          .VERSION KERNSUN  1.00  881114
* ORIG. 14/11/88, JZ
*/
#include <setjmp.h>
qnexte_()
{     static jmp_buf  myenv;
      static int ireent = 0;
      static int j = 7;

      if (ireent)  longjmp(myenv, j);

      ireent = 77;
      setjmp(myenv);
      qnext_();
}
