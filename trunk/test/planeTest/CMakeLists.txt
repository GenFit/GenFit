cmake_minimum_required(VERSION 2.6)

IF(NOT DEFINED ENV{ROOTSYS})
   MESSAGE(FATAL_ERROR "You did not set the environment variable ROOTSYS. Please check your ROOT configuration.") 
ENDIF(NOT DEFINED ENV{ROOTSYS})

IF(NOT DEFINED ENV{GENFIT})
   MESSAGE(FATAL_ERROR "You did not set the environment variable GENFIT.") 
ENDIF(NOT DEFINED ENV{GENFIT})

SET(CMAKE_CXX_FLAGS "-g -O2")
#SET(CMAKE_VERBOSE_MAKEFILE ON)

Execute_process(COMMAND root-config --glibs --cflags
	OUTPUT_VARIABLE ROOT_LIBS
    	OUTPUT_STRIP_TRAILING_WHITESPACE)

INCLUDE_DIRECTORIES($ENV{ROOTSYS}/include $ENV{GENFIT}/core)

ADD_EXECUTABLE(main main.cxx)
TARGET_LINK_LIBRARIES(main ${ROOT_LIBS} -L$ENV{GENFIT}/lib -lGeom -lgenfit)

ADD_EXECUTABLE(mainRead mainRead.cxx)
TARGET_LINK_LIBRARIES(mainRead ${ROOT_LIBS} -L$ENV{GENFIT}/lib -lGeom -lgenfit)


