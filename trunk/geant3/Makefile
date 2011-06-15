# $Id: Makefile 247 2009-10-09 13:20:34Z brun $

############################### geant321 Makefile #############################

PACKAGE   = geant321

ifeq ($(PLATFORM),)
PLATFORM = $(shell root-config --arch)
endif

TOPDIR  = $(shell pwd)
BINDIR  = $(TOPDIR)/tgt_$(PLATFORM)
LIBDIR  = $(TOPDIR)/lib/tgt_$(PLATFORM)

ifeq ($(ROOTSYS),)
ROOT_INCDIR = $(shell root-config --incdir)
ROOT_BINDIR = $(shell root-config --prefix)/bin
ROOT_ETCDIR = $(shell root-config --etcdir)/vmc
else
ROOT_INCDIR = $(ROOTSYS)/include
ROOT_BINDIR = $(ROOTSYS)/bin
ROOT_ETCDIR = $(ROOTSYS)/etc/vmc
endif

include $(ROOT_ETCDIR)/Makefile.$(PLATFORM)

ifneq ($(findstring gcc,$(CC)),)
GCC_MAJOR    := $(shell $(CC) -dumpversion 2>&1 | cut -d'.' -f1)
GCC_MINOR    := $(shell $(CC) -dumpversion 2>&1 | cut -d'.' -f2)
GCC_VERS     := gcc-$(GCC_MAJOR).$(GCC_MINOR)
endif


############################### Sources #######################################

GDIRS:=	added gbase gcons geocad ggeom gheisha ghits ghrout ghutils \
	giface giopa gkine gparal gphys gscan gstrag gtrak matx55 miface \
	miguti neutron peanut fiface cgpack fluka block comad erdecks erpremc \
        minicern gdraw

include $(ROOT_ETCDIR)/MakeRules


# C++ Headers

DHDRS          = TGeant3/TGeant3.h 


# Dummy Geant Objects

DOBJS          = $(patsubst %.cxx,$(BINDIR)/%.o,$(DSRCS))

# Library dictionary

GDICT    := $(BINDIR)/TGeant3/geant3Cint.cxx
GDICTH   := $(DICT:.cxx=.h)
GDICTO   := $(patsubst %.cxx,%.o,$(GDICT))

# Sources

FSRC	:= $(wildcard $(patsubst %,%/*.F,$(GDIRS))) gcinit.F
FSRC	:= $(filter-out gtrak/grndm%.F,$(FSRC))
FSRC	:= $(filter-out erdecks/eustep.F,$(FSRC))

ifeq ($(PLATFORM),linux)
	  FSRC += minicern/lnxgs/rdmin.F
endif
ifneq ($(PLATFORM),solarisCC5)
	  FSRC := $(filter-out minicern/uset.F,$(FSRC))
endif
CSRC	:= $(wildcard $(patsubst %,%/*.c,$(GDIRS))) 
CSRC     := $(filter-out minicern/lnblnk.c,$(CSRC)) 
ifeq ($(PLATFORM),linux)
	  CSRC += minicern/lnxgs/ishftr.c
endif
ifeq ($(PLATFORM),macosx)
	  CSRC += minicern/lnxgs/ishftr.c
endif
ifeq ($(PLATFORM),macosx64)
	  CSRC += minicern/lnxgs/ishftr.c
endif
ifeq ($(PLATFORM),linuxicc)
	  CSRC += minicern/lnxgs/ishftr.c
endif
ifeq ($(PLATFORM),linuxia64ecc)
	  CSRC += minicern/lnxgs/ishftr.c
endif
ifeq ($(PLATFORM),linuxx8664gcc)
	  CSRC += minicern/lnxgs/ishftr.c
endif
ifeq ($(PLATFORM),linuxia64gcc)
          CSRC += minicern/lnxgs/ishftr.c
endif
ifeq ($(PLATFORM),hpuxacc)
	  CSRC += minicern/hpxgs/traceqc.c
endif
ifeq ($(PLATFORM),hpuxacc)
	  CSRC += minicern/lnblnk.c
endif
ifeq ($(PLATFORM),macosx)
	  CSRC += minicern/lnblnk.c
endif
ifeq ($(PLATFORM),macosx64)
	  CSRC += minicern/lnblnk.c
endif
CSRC	:= $(filter-out ,$(CSRC))

CXXSRC	:= $(wildcard $(patsubst %,%/*.cxx,$(GDIRS))) \
           $(wildcard TGeant3/*.cxx)
SRCS	:= $(FSRC) $(CSRC) $(CXXSRC)

# C++ Headers

HDRS    := $(CXXSRC:.cxx=.h) TGeant3/geant3LinkDef.h
HDRS	:= $(filter-out comad/gcadd.h,$(HDRS))  

# Objects

FOBJ	:= $(patsubst %.F,$(BINDIR)/%.o,$(FSRC))
COBJ	:= $(patsubst %.c,$(BINDIR)/%.o,$(CSRC))
CXXOBJ	:= $(patsubst %.cxx,$(BINDIR)/%.o,$(CXXSRC))
OBJS	:= $(FOBJ) $(COBJ) $(CXXOBJ) $(GDICTO)


# C++ compilation flags

CXXFLAGS := $(CXXOPTS) $(CLIBCXXOPTS) -I. -I$(ROOT_INCDIR) -ITGeant3 -fPIC

# C compilation flags

CFLAGS      := $(COPT) $(CLIBCOPT) -I. -Iminicern -fPIC

# FORTRAN compilation flags

FFLAGS      := $(FOPT) $(CLIBFOPT) -I. -Iminicern -fPIC
# ifeq ($(PLATFORM),linux)
#    FFLAGS      := $(filter-out -O%,$(FFLAGS))  
# endif

DEPINC 		+= -I. -I$(ROOT_INCDIR)

############################### Targets #######################################


SLIBRARY	= $(LIBDIR)/lib$(PACKAGE).$(SL)
ALIBRARY	= $(LIBDIR)/lib$(PACKAGE).a

ifeq ($(PLATFORM),OSF1)
        default:	depend $(ALIBRARY) $(SLIBRARY) 
else
        default:	depend $(SLIBRARY)
endif

$(LIBDIR)/lib$(PACKAGE).$(SL):  $(OBJS)
$(LIBDIR)/lib$(PACKAGE).a:  $(OBJS)

DICT:= $(GDICT) $(DDICT)

$(GDICT): $(HDRS)

$(DDICT): $(DHDRS)

depend:		$(SRCS)

TOCLEAN		= $(BINDIR)
TOCLEANALL	= $(BINDIR) $(LIBDIR)

MAKEDIST	= $(ROOT_ETCDIR)/g3_makedist.sh $(GCC_VERS) lib
MAKEDISTSRC	= $(ROOT_ETCDIR)/g3_makedist.sh $(GCC_VERS)

include $(ROOT_ETCDIR)/MakeMacros

############################### Dependencies ##################################

-include $(BINDIR)/Make-depend
