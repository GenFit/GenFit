cmake_minimum_required(VERSION 2.6)


#check for ENV variables that are needed
IF(NOT DEFINED ENV{ROOTSYS})
   MESSAGE(FATAL_ERROR "You did not set the environment variable ROOTSYS. Please check your ROOT configuration.") 
ENDIF(NOT DEFINED ENV{ROOTSYS})
#IF(NOT DEFINED ENV{GENFIT})
#   MESSAGE(FATAL_ERROR "You did not set the environment variable GENFIT.") 
#ENDIF(NOT DEFINED ENV{GENFIT})
IF(NOT DEFINED ENV{G4INCLUDE})
   MESSAGE(FATAL_ERROR "You did not set the environment variable G4INCLUDE. Please check your Geant4 configuration.") 
ENDIF(NOT DEFINED ENV{G4INCLUDE})

INCLUDE_DIRECTORIES(
  $ENV{ROOTSYS}/include
  ${GENFIT}/core
  ${GENFIT}/G4eTrackRep
  $ENV{G4INCLUDE}
  $ENV{G4INCLUDE}/geant4
  $ENV{CLHEP_INCLUDE_DIR}
)

FILE(GLOB _ALL_FILES $ENV{CLHEP_INCLUDE_DIR} $ENV{CLHEP_INCLUDE_DIR}/*)
FOREACH(_FILE ${_ALL_FILES})
  IF(IS_DIRECTORY ${_FILE})
    INCLUDE_DIRECTORIES(${_FILE})
  ENDIF()
ENDFOREACH()


AUX_SOURCE_DIRECTORY(
  $ENV{GENFIT}/G4eTrackRep 
  GENFIT_G4e_SRCS
)

SET(GENFIT_G4e_HEADERS 
  G4eTrackRep.h
)

set(GENFIT_G4e_LINKDEF  genfitG4eLinkDef.h)
set(GENFIT_G4e_DICTIONARY genfitG4eDict.cxx) 

FIND_PROGRAM(ROOT_CINT_EXECUTABLE
  NAMES rootcint
  PATHS $ENV{ROOTSYS}/bin
  NO_DEFAULT_PATH)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_BINARY_DIR}/${GENFIT_G4e_DICTIONARY}
  COMMAND ${ROOT_CINT_EXECUTABLE}
    ARGS -f ${CMAKE_BINARY_DIR}/${GENFIT_G4e_DICTIONARY} -c
      -DHAVE_CONFIG_H
      -I${GENFIT}/core
      -I${GENFIT}/G4eTrackRep
      -I$ENV{G4INCLUDE}
      -I$ENV{G4INCLUDE}/geant4
      -I$ENV{CLHEP_INCLUDE_DIR}
      ${GENFIT_G4e_HEADERS}
      ${GENFIT_G4e_LINKDEF}
    DEPENDS
      ${GENFIT_G4e_HEADERS}
      ${GENFIT_G4e_LINKDEF}
)

SET(GENFIT_G4e_SRCS ${GENFIT_G4e_SRCS} ${CMAKE_BINARY_DIR}/${GENFIT_G4e_DICTIONARY})


################ install ###################

#add ROOT lib
Execute_process(COMMAND root-config --glibs --cflags
	OUTPUT_VARIABLE ROOT_LIBS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)



ADD_LIBRARY(genfitG4e SHARED ${GENFIT_G4e_SRCS})
TARGET_LINK_LIBRARIES(
  genfitG4e 
  ${ROOT_LIBS} 
  -L${GENFITLIBDIR} 
  -lgenfit
  -lGeom 
  -lVMC 
  -lEG
  -L$ENV{G4LIB}
  -L$ENV{G4LIB}/$ENV{G4SYSTEM}
  -lG4FR             -lG4intercoms   -lG4particles     -lG4readout   -lG4visHepRep
  -lG3toG4               -lG4geometry       -lG4interfaces  -lG4persistency   -lG4run       -lG4vis_management
  -lG4digits_hits        -lG4global         -lG4materials   -lG4physicslists  -lG4tracking  -lG4visXXX
  -lG4error_propagation  -lG4GMocren        -lG4modeling    -lG4processes     -lG4track     -lG4VRML
  -lG4event              -lG4graphics_reps  -lG4parmodels   -lG4RayTracer     -lG4Tree      -lG4zlib
  -L$ENV{CLHEP_LIB_DIR}
  -lCLHEP-Evaluator-2.1.0.1     -lCLHEP-Geometry-2.1.0.1           -lCLHEP-RandomObjects-2.1.0.1
  -lCLHEP-2.1.0.1               -lCLHEP-Exceptions-2.1.0.1         -lCLHEP-Matrix-2.1.0.1           -lCLHEP-RefCount-2.1.0.1
  -lCLHEP-Cast-2.1.0.1          -lCLHEP-GenericFunctions-2.1.0.1   -lCLHEP-Random-2.1.0.1           -lCLHEP-Vector-2.1.0.1
)
INSTALL(TARGETS genfitG4e DESTINATION ${GENFITLIBINSTALLDIR} )

