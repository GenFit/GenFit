#if 0
* This pilot patch was created from kerngen.car patch _kerngen
* This directory was created from kernfor.car patch qdefault
*   Character set is ASCII
*  Internal double-precision
*            copy vectors as floating normally
#endif
#if defined(CERNLIB_LXIA64)
#if 0
*
* Take normal LINUX as basis for Itanium
#endif
#ifndef CERNLIB_LINUX
#define CERNLIB_LINUX
#endif
#endif
#if defined(CERNLIB_MSDOS)
#include "kerngen/kdos.h"
#endif
#if defined(CERNLIB_WINNT)
#include "kerngen/kwnt.h"
#endif
#if (defined(CERNLIB_DECS))&&(!defined(CERNLIB_QMVAOS))
#include "kerngen/kvmi.h"
#endif
#if (defined(CERNLIB_DECS))&&(defined(CERNLIB_QMVAOS))
#include "kerngen/kvaos.h"
#endif
#if defined(CERNLIB_HPUX)
#include "kerngen/khpx.h"
#endif
#if defined(CERNLIB_IBMRT)
#include "kerngen/kirt.h"
#endif
#if defined(CERNLIB_LINUX)
#include "kerngen/klnx.h"
#endif
#if defined(CERNLIB_MACMPW)
#include "kerngen/kmpw.h"
#endif
#if defined(CERNLIB_OS9)
#if 0
* Added at release 94B
#endif
#include "kerngen/kos9.h"
#endif
#if defined(CERNLIB_SGI)
#include "kerngen/ksgi.h"
#endif
#if defined(CERNLIB_SUN)
#include "kerngen/ksun.h"
#endif
#if (defined(CERNLIB_VAXVMS))&&(!defined(CERNLIB_QMALPH))
#include "kerngen/kvax.h"
#endif
#if (defined(CERNLIB_VAXVMS))&&(defined(CERNLIB_QMALPH))
#include "kerngen/kalph.h"
#endif
#ifndef CERNLIB_A4
#define CERNLIB_A4
#endif
#ifndef CERNLIB_B32
#define CERNLIB_B32
#endif
#ifndef CERNLIB_HEX
#define CERNLIB_HEX
#endif
#if !defined(CERNLIB_QEBCDIC)
#ifndef CERNLIB_QASCII
#define CERNLIB_QASCII
#endif
#endif
#if defined(CERNLIB_B32)||defined(CERNLIB_B36)
#ifndef CERNLIB_INTDOUBL
#define CERNLIB_INTDOUBL
#endif
#endif
#if defined(CERNLIB_QX_SC)
#ifdef CERNLIB_QXNO_SC
#undef CERNLIB_QXNO_SC
#endif
#endif
#if defined(CERNLIB_QXNO_SC)
#ifdef CERNLIB_QX_SC
#undef CERNLIB_QX_SC
#endif
#endif
#if defined(CERNLIB_QXNO_SC)||defined(CERNLIB_QX_SC)
#ifdef CERNLIB_QXCAPT
#undef CERNLIB_QXCAPT
#endif
#endif
#if (!defined(CERNLIB_QXNO_SC))&&(!defined(CERNLIB_QXCAPT))
#ifndef CERNLIB_QX_SC
#define CERNLIB_QX_SC
#endif
#endif
#if defined(CERNLIB_SHIFT)
#ifndef CERNLIB_PROJSHIFT
#define CERNLIB_PROJSHIFT
#endif
#endif

#ifndef type_of_call
#define type_of_call
#endif

