cmake_minimum_required(VERSION 2.6)

IF(NOT DEFINED ENV{ROOTSYS})
   MESSAGE(FATAL_ERROR "You did not set the environment variable ROOTSYS. Please check your ROOT configuration.") 
ENDIF(NOT DEFINED ENV{ROOTSYS})

SET(CMAKE_CXX_FLAGS "-g -O2 -Wall -Wshadow -ftree-vectorize -ffast-math -march=native -ftree-vectorizer-verbose=0")

include(macros.cmake)

#SET(CMAKE_VERBOSE_MAKEFILE ON)



SET(GENFIT $ENV{GENFIT})
SET(VMC $ENV{VMC})

SET(GENFITLIBDIR ${GENFIT}/lib)
SET(GENFITLIBINSTALLDIR ${GENFITLIBDIR})
SET(LIBRARY_OUTPUT_PATH ${GENFITLIBDIR})


add_subdirectory(core)
add_subdirectory(GeaneTrackRep2)
add_subdirectory(RKTrackRep)
add_subdirectory(GFFastKalman)
#add_subdirectory(RKTrackRepXY)
#add_subdirectory(LSLtrackRep)
#add_subdirectory(SlTrackRep)
#add_subdirectory(G4eTrackRep)
#add_subdirectory(GFRave)
#add_subdirectory(tool)

