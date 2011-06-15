#if defined(CERNLIB_LXIA64)
*
* Take normal LINUX as basis for Itanium
#ifndef CERNLIB_LINUX
#define CERNLIB_LINUX
#endif
#endif

#if defined(CERNLIB_LINUX) || defined(CERNLIB_HPUX) || defined(CERNLIB_SUN)
#if !defined(CERNLIB_UNIX)
#define CERNLIB_UNIX
#endif
#endif

#if defined(CERNLIB_UNIX)
#define CERNLIB_Z32
#define CERNLIB_QMUIX
#define CERNLIB_A4
#define CERNLIB_B32
#define CERNLIB_HEX
#define CERNLIB_QHOLL
# ifndef CERNLIB_QFEPC
#define CERNLIB_EQUHOLCH
# endif
#define CERNLIB_QTRHOLL
#define CERNLIB_QASCII
#define CERNLIB_QCFIO
#define CERNLIB_QIEEE
#define CERNLIB_FQXISN
#define CERNLIB_QPRINT
#define CERNLIB_QDEBUG
#define CERNLIB_QDEBPRI
#define CERNLIB_QISASTD
#define CERNLIB_QMILSTD
#define CERNLIB_FZALFA
#define CERNLIB_FZCHANNEL
#define CERNLIB_FZDACC
#define CERNLIB_FZDACCF
#define CERNLIB_FZDACCH
#define CERNLIB_FZDACCL
#define CERNLIB_FZFFNAT
#define CERNLIB_FZFORTRAN
#define CERNLIB_FZLIBC
#define CERNLIB_FZMEMORY
#define CERNLIB_JZTIME
#define CERNLIB_RZBYTES
#define CERNLIB_RZFRECL
#endif

#if defined(CERNLIB_IBMRT)
#define CERNLIB_QMIRT 
#undef CERNLIB_QMILSTD
#endif

#if defined(CERNLIB_HPUX)
#define CERNLIB_QMHPX
#endif

#if defined(CERNLIB_SUN)
#define CERNLIB_QMSUN
#if 0
* CERNLIB_BUGLRSHFT to get round the lrshft bug in Sun f77 3.0.x
* At CERN a patch was applied for Solaris only
*
#endif
#ifndef CERNLIB_SOLARIS
#define CERNLIB_BUGLRSHFT
#endif

#undef CERNLIB_QMILSTD
#undef CERNLIB_QISASTD
#endif

#if defined(CERNLIB_SGI)
#define CERNLIB_QMSGI
#undef CERNLIB_RZBYTES
#endif

#if (defined(CERNLIB_DECS))||(defined(CERNLIB_DECOSF))
#define CERNLIB_QMVMI
#define CERNLIB_FQNEEDCV
#undef CERNLIB_Z32
#undef CERNLIB_FQXISN
#undef CERNLIB_RZBYTES
#undef CERNLIB_QMILSTD
#endif

#if defined(CERNLIB_WINNT)
#define CERNLIB_QMDOS
#define CERNLIB_QS_WNT
#ifndef _X86_
# define CERNLIB_QFDEC
#endif
#define CERNLIB_FQNEEDCV
#define CERNLIB_F77TRARG
#undef CERNLIB_Z32
#undef CERNLIB_FQXISN
#ifdef CERNLIB_QFDEC
# undef CERNLIB_RZBYTES
#endif
#endif

#if defined(CERNLIB_LINUX)
#define CERNLIB_QMLNX
#if !defined(CERNLIB_PPC)
# define CERNLIB_FQNEEDCV
# undef CERNLIB_Z32
# undef CERNLIB_FQXISN
#endif
#endif

#if defined(CERNLIB_VAXVMS)
#define CERNLIB_QMVAX
#define CERNLIB_A4
#define CERNLIB_B32
#define CERNLIB_HEX
#define CERNLIB_QHOLL
#define CERNLIB_EQUHOLCH
#define CERNLIB_QTRHOLL
#define CERNLIB_QASCII
#define CERNLIB_QCFIO
#define CERNLIB_QMILSTD
#define CERNLIB_QPRINT
#define CERNLIB_QDEBUG
#define CERNLIB_QDEBPRI
#define CERNLIB_FZALFA
#define CERNLIB_FZCHANNEL
#define CERNLIB_FZDACC
#define CERNLIB_FZDACCF
#define CERNLIB_FZDACCH
#define CERNLIB_FZDACCL
#define CERNLIB_FZFFNAT
#define CERNLIB_FZFORTRAN
#define CERNLIB_FZLIBC
#define CERNLIB_FZMEMORY
#define CERNLIB_JZTIME
#define CERNLIB_RZFRECL
#define CERNLIB_FQIE3FDC
#define CERNLIB_FQIE3FSC
#define CERNLIB_FQIE3TDC
#define CERNLIB_FQIE3TSC
#define CERNLIB_FQNEEDCV
#endif
