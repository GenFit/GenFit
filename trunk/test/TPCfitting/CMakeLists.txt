cmake_minimum_required(VERSION 2.6)

IF(NOT DEFINED ENV{ROOTSYS})
   MESSAGE(FATAL_ERROR "You did not set the environment variable ROOTSYS. Please check your ROOT configuration.") 
ENDIF(NOT DEFINED ENV{ROOTSYS})

IF(NOT DEFINED ENV{GENFIT})
   MESSAGE(FATAL_ERROR "You did not set the environment variable GENFIT.") 
ENDIF(NOT DEFINED ENV{GENFIT})

IF(NOT DEFINED ENV{VMC})
   MESSAGE(FATAL_ERROR "You did not set the environment variable VMC for your virtual Monte Carlo installation.") 
ENDIF(NOT DEFINED ENV{VMC})


SET(CMAKE_CXX_FLAGS "-g -O2")
#SET(CMAKE_VERBOSE_MAKEFILE ON)

Execute_process(COMMAND root-config --glibs --cflags
	OUTPUT_VARIABLE ROOT_LIBS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

INCLUDE_DIRECTORIES($ENV{ROOTSYS}/include $ENV{GENFIT}/core
$ENV{GENFIT}/GeaneTrackRep2 $ENV{GENFIT}/RKTrackRep $ENV{VMC}/TGeant3)

SET(MYDICT myDict.cxx)
SET(MYLINKDEF linkdef.h)
SET(DICTHEADERS StripHit.h PointHit.h)
FIND_PROGRAM(ROOT_CINT_EXECUTABLE
	NAMES rootcint
	PATHS $ENV{ROOTSYS}/bin
	NO_DEFAULT_PATH)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_BINARY_DIR}/${MYDICT}
	COMMAND ${ROOT_CINT_EXECUTABLE}
      	ARGS -f ${CMAKE_BINARY_DIR}/${MYDICT} -c -I$ENV{GENFIT}/core -I$ENV{GENFIT}/GeaneTrackRep  -DHAVE_CONFIG_H
      	${DICTHEADERS} ${MYLINKDEF}
      	DEPENDS ${DICTHEADERS} ${MYLINKDEF})

ADD_EXECUTABLE(main main.cxx StripHit.cxx PointHit.cxx ${MYDICT})
TARGET_LINK_LIBRARIES(main ${ROOT_LIBS} -lGeom -L$ENV{GENFIT}/lib
-lgenfit -lgenfitGeane -lgenfitRK)

