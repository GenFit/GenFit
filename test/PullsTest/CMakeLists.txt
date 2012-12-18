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

Execute_process(COMMAND root-config --glibs --cflags --evelibs
	OUTPUT_VARIABLE ROOT_LIBS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

Execute_process(COMMAND root-config --cflags
	OUTPUT_VARIABLE ROOT_CFLAGS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

SET(CMAKE_CXX_FLAGS "-g -O2 ${ROOT_CFLAGS}")
#SET(CMAKE_VERBOSE_MAKEFILE ON)

INCLUDE_DIRECTORIES(
	$ENV{GENFIT}/core 
	$ENV{GENFIT}/RecoHitExamples 
	$ENV{GENFIT}/RKTrackRep 
	$ENV{GENFIT}/GeaneTrackRep2 
	$ENV{VMC}/TGeant3 
	$ENV{GENFIT}/eventDisplay
)

SET(MYDICT myDict.cxx)
SET(MYLINKDEF linkdef.h)
SET(DICTHEADERS 
	$ENV{GENFIT}/eventDisplay/GenfitDisplay.h
)
FIND_PROGRAM(ROOT_CINT_EXECUTABLE
	NAMES rootcint
	PATHS $ENV{ROOTSYS}/bin
	NO_DEFAULT_PATH)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_BINARY_DIR}/${MYDICT}
	COMMAND ${ROOT_CINT_EXECUTABLE}
      	ARGS -f ${CMAKE_BINARY_DIR}/${MYDICT} -c -I$ENV{GENFIT}/core -I$ENV{GENFIT}/RecoHitExamples -I$ENV{GENFIT}/GeaneTrackRep2 -I${VMC}/TGeant3 -DHAVE_CONFIG_H
      	${DICTHEADERS} ${MYLINKDEF}
      	DEPENDS ${DICTHEADERS} ${MYLINKDEF})

ADD_EXECUTABLE(main 
	main.cxx 
	$ENV{GENFIT}/eventDisplay/GenfitDisplay.h 
	$ENV{GENFIT}/eventDisplay/GenfitDisplay.cxx 
	${MYDICT}
)
TARGET_LINK_LIBRARIES(main 
	${ROOT_LIBS} 
	-lGeom 
	-lVMC 
	-lEG
	-L$ENV{GENFIT}/lib
	-lgenfit 
	-lgenfitRK 
	-lgenfitGeane
	-lmyRecoHitExamples
)

