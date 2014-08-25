# #require a certain version
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

# declare the project name
PROJECT(genfit2)

# set project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 2 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

# install destinations can be passed via the command line:
# cmake -DLIBRARY_OUTPUT_DIRECTORY:PATH=.... -DINCLUDE_OUTPUT_DIRECTORY:PATH=....
# otherwise, set default install destinations
IF(NOT DEFINED LIBRARY_OUTPUT_DIRECTORY)
	set( LIBRARY_OUTPUT_DIRECTORY ./lib )
ENDIF()
IF(NOT DEFINED INCLUDE_OUTPUT_DIRECTORY)
	set( INCLUDE_OUTPUT_DIRECTORY ./include )
ENDIF()

# define path to macros and settings
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
# include the macros and settings
INCLUDE(settings)
INCLUDE(macros)

FIND_PACKAGE(ROOT REQUIRED COMPONENTS Geom)

# either, the environment variable RAVEPATH has to be specified, and RAVE's pkgConfig will be used to determine everything we need
# or, RAVE_LDFLAGS, RAVE_INCLUDE_DIRS and RAVE_CFLAGS have to be set via the command-line


IF(DEFINED ENV{RAVEPATH})
	MESSAGE("Environment variable RAVEPATH is set. Try to build GFRave.") 
	# get compiler flags from rave
	INCLUDE(FindPkgConfig)
	SET(ENV{PKG_CONFIG_PATH} $ENV{RAVEPATH})
	pkg_search_module(RAVE rave REQUIRED)
	
	SET(RAVE True)
ELSEIF(DEFINED RAVE_LDFLAGS)
	IF(DEFINED RAVE_INCLUDE_DIRS)
		IF(DEFINED RAVE_CFLAGS)
			MESSAGE("RAVE_LDFLAGS, RAVE_INCLUDE_DIRS and RAVE_CFLAGS are set. Try to build GFRave.") 
			SET(RAVE True)
		ENDIF()
	ENDIF()
ELSE()
   	MESSAGE("You did not define the environment variable RAVEPATH which is needed to find rave. GFRave will not be built.") 
ENDIF()

IF(DEFINED RAVE)
	# turn comma separated lists into space separated strings
	string (REPLACE ";" " " RAVE_LDFLAGS_STR "${RAVE_LDFLAGS}")
	string (REPLACE ";" " " RAVE_INCLUDE_DIRS_STR "${RAVE_INCLUDE_DIRS}")
	string (REPLACE ";" " " RAVE_CFLAGS_STR "${RAVE_CFLAGS}")

	SET (CMAKE_CXX_FLAGS ${RAVE_CFLAGS_STR})
	
	SET(GF_INC_DIRS
		./GFRave/include/
	)
	
	AUX_SOURCE_DIRECTORY( ./GFRave/src  library_sources )
ENDIF()

# Set GENFIT include directories
SET(GF_INC_DIRS
	./core/include/ 
	./eventDisplay/include/
	./fields/include/
	./finitePlanes/include/ 
	./fitters/include/
	./GBL/include/
	./measurements/include/ 
	./trackReps/include/
	./utilities/include/
    ${GF_INC_DIRS}
)

#search for boost, provide some information to user
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    MESSAGE("-- Boost is found at: ${Boost_INCLUDE_DIR}")
    INSTALL( DIRECTORY ${GF_INC_DIRS}
    	DESTINATION ${INCLUDE_OUTPUT_DIRECTORY}
	    PATTERN ".svn" EXCLUDE 
	 )
ELSE()
    MESSAGE( " Boost was not found on your system, using the one shipped with GENFIT.")
    INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/extern/ )
    INSTALL( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/extern/ ${GF_INC_DIRS}  
    	DESTINATION ${INCLUDE_OUTPUT_DIRECTORY}
    	PATTERN ".svn" EXCLUDE 
    )
ENDIF()


IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Uncomment for coverage analysis.
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")

include(CheckCXXCompilerFlag)
SET(CPP_STANDARD c++03)
# ROOT 6 requires C++11 support, enforced via headers
IF (ROOT_found_version GREATER 59900)
    SET(CPP_STANDARD c++11)
ENDIF()
MESSAGE(STATUS "Using C++ standard ${CPP_STANDARD}.")
check_cxx_compiler_flag(-std=${CPP_STANDARD} HAS_STD_FLAG)
check_cxx_compiler_flag(--std=${CPP_STANDARD} HAS_STD_FLAG_ALTERNATIVE)
if (HAS_STD_FLAG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=${CPP_STANDARD}")
elseif (HAS_STD_FLAG_ALTERNATIVE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=${CPP_STANDARD}")
else()
  message("Neither -std=${CPP_STANDARD} nor --std=${CPP_STANDARD} seem to work.  Moving on.")
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wno-long-long -Wshadow -Werror=overloaded-virtual"   )
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g3")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O3 -g3")




# include directories 
INCLUDE_DIRECTORIES( 
	BEFORE 
	${GF_INC_DIRS}
	${ROOT_INCLUDE_DIRS})


# add a target to generate API documentation with Doxygen
FIND_PACKAGE(Doxygen)

if(DOXYGEN_FOUND)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)

    ADD_CUSTOM_TARGET(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM)
ELSE()
    MESSAGE( STATUS "Doxygen not installed or found -- API documentation can't be generated!" )
ENDIF(DOXYGEN_FOUND)


# declare the ROOT library path
LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})


# add the list of source code files -- careful, when a new file is added!
AUX_SOURCE_DIRECTORY( ./core/src          library_sources )
AUX_SOURCE_DIRECTORY( ./eventDisplay/src  library_sources )
AUX_SOURCE_DIRECTORY( ./fields/src        library_sources )
AUX_SOURCE_DIRECTORY( ./finitePlanes/src  library_sources )
AUX_SOURCE_DIRECTORY( ./fitters/src       library_sources )
AUX_SOURCE_DIRECTORY( ./GBL/src           library_sources )
AUX_SOURCE_DIRECTORY( ./measurements/src  library_sources )
AUX_SOURCE_DIRECTORY( ./trackReps/src     library_sources )
AUX_SOURCE_DIRECTORY( ./utilities/src     library_sources )


# Dictionary generation.  For the time being, we list classes one-by-one.
MACRO (ROOT_GENERATE_DICTIONARY LIBNAME INFILES INCLUDE_DIRS_IN LINKDEF_FILE OUTFILE)
  set (INCLUDE_DIRS)
  foreach (_current_FILE ${INCLUDE_DIRS_IN})
    set(INCLUDE_DIRS ${INCLUDE_DIRS} -I${_current_FILE})   
  endforeach (_current_FILE ${INCLUDE_DIRS_IN})

  SET(EXTRA_DICT_ARGS "")
  IF(ROOT_found_version GREATER 59900)
    STRING(REGEX REPLACE "^(.*)\\.(.*)$" "\\1_rdict.pcm" OUTFILE_PCM "${OUTFILE}")
    STRING(REGEX REPLACE "^(.*)Dict\\.(.*)$" "\\1.rootmap" OUTFILE_RMF "${OUTFILE}")
    SET (OUTFILES ${OUTFILE} ${OUTFILE_PCM} ${OUTFILE_RMF})
    SET(EXTRA_DICT_ARGS -inlineInputHeader -rmf ${OUTFILE_RMF} -rml ${LIBNAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
  ELSE()
    STRING(REGEX REPLACE "^(.*)\\.(.*)$" "\\1.h" OUTFILE_H "${OUTFILE}")
    SET (OUTFILES ${OUTFILE} ${OUTFILE_H})
  ENDIF()
  
  ADD_CUSTOM_COMMAND(OUTPUT ${OUTFILES}
    COMMAND ${ROOT_CINT_EXECUTABLE}
    ARGS -f ${OUTFILE} ${EXTRA_DICT_ARGS} -c ${INCLUDE_DIRS} ${INFILES} ${LINKDEF_FILE}
    DEPENDS ${INFILES} ${LINKDEF_FILE})

  STRING(REGEX REPLACE "/" "" OUTFILE_NO_SLASH "${OUTFILE}")
  ADD_CUSTOM_TARGET(generate_${OUTFILE_NO_SLASH} DEPENDS ${OUTFILE})

  IF(ROOT_found_version GREATER 59900)
    # PCM files and rootmap-files should be installed
    INSTALL(FILES ${OUTFILE_RMF} ${OUTFILE_PCM}
            DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
  ENDIF()
ENDMACRO (ROOT_GENERATE_DICTIONARY)

SET(CORE_DICTIONARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsFinitePlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsFitterInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsHMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsTrackRep.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/DetPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/FitStatus.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/MaterialProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/MeasuredStateOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/MeasurementFactory.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/MeasurementOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/MeasurementProducer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/StateOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/ThinScatterer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/Track.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/TrackCand.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/TrackCandHit.h
  ${CMAKE_CURRENT_SOURCE_DIR}/core/include/TrackPoint.h
  ${CMAKE_CURRENT_SOURCE_DIR}/finitePlanes/include/RectangularFinitePlane.h
)
ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
                         "${CORE_DICTIONARY_SOURCES}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/core/include"
			 "${CMAKE_CURRENT_SOURCE_DIR}/core/src/coreLinkDef.h"
			 "${CMAKE_CURRENT_BINARY_DIR}/coreRootDict.cc")
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/coreRootDict.cc" "${library_sources}")

SET(EVENTDISPLAY_DICTIONARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/eventDisplay/include/EventDisplay.h
)
ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
                         "${EVENTDISPLAY_DICTIONARY_SOURCES}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/fitters/include;${CMAKE_CURRENT_SOURCE_DIR}/eventDisplay/include"
			 "${CMAKE_CURRENT_SOURCE_DIR}/eventDisplay/src/eventDisplayLinkDef.h"
			 "${CMAKE_CURRENT_BINARY_DIR}/eventDisplayRootDict.cc")
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/eventDisplayRootDict.cc" "${library_sources}")

SET(FITTERS_DICTIONARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/AbsKalmanFitter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFittedStateOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFitStatus.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/ReferenceStateOnPlane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFitterRefTrack.h

  ${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFitter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/KalmanFitterInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fitters/include/DAF.h
)
ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
                         "${FITTERS_DICTIONARY_SOURCES}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/fitters/include"
			 "${CMAKE_CURRENT_SOURCE_DIR}/fitters/src/fittersLinkDef.h"
			 "${CMAKE_CURRENT_BINARY_DIR}/fittersRootDict.cc")
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/fittersRootDict.cc" "${library_sources}")

SET(GBL_DICTIONARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GFGbl.h
  ${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblFitter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/ICalibrationParametersDerivatives.h
  ${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblFitterInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblFitStatus.h
)
ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
                         "${GBL_DICTIONARY_SOURCES}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/GBL/include"
			 "${CMAKE_CURRENT_SOURCE_DIR}/GBL/src/GBLLinkDef.h"
			 "${CMAKE_CURRENT_BINARY_DIR}/GBLRootDict.cc")
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/GBLRootDict.cc" "${library_sources}")

SET(MEASUREMENTS_DICTIONARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixPhi.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixU.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixUnit.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixV.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/HMatrixUV.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/FullMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/PlanarMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/ProlateSpacepointMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/SpacepointMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WireMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WirePointMeasurement.h
  ${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WireTrackCandHit.h
)
ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
                         "${MEASUREMENTS_DICTIONARY_SOURCES}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/measurements/include"
			 "${CMAKE_CURRENT_SOURCE_DIR}/measurements/src/measurementsLinkDef.h"
			 "${CMAKE_CURRENT_BINARY_DIR}/measurementsRootDict.cc")
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/measurementsRootDict.cc" "${library_sources}")

SET(TRACKREPS_DICTIONARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/trackReps/include/RKTrackRep.h
)
ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
                         "${TRACKREPS_DICTIONARY_SOURCES}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/trackReps/include"
			 "${CMAKE_CURRENT_SOURCE_DIR}/trackReps/src/trackRepsLinkDef.h"
			 "${CMAKE_CURRENT_BINARY_DIR}/trackRepsRootDict.cc")
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/trackRepsRootDict.cc" "${library_sources}")

SET(UTILITIES_DICTIONARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/HelixTrackModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/MeasurementCreator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/mySpacepointDetectorHit.h
  ${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/mySpacepointMeasurement.h
)
ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
                         "${UTILITIES_DICTIONARY_SOURCES}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/utilities/include;${CMAKE_CURRENT_SOURCE_DIR}/measurements/include"
			 "${CMAKE_CURRENT_SOURCE_DIR}/utilities/src/utilitiesLinkDef.h"
			 "${CMAKE_CURRENT_BINARY_DIR}/utilitiesRootDict.cc")
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/utilitiesRootDict.cc" "${library_sources}")

if(DEFINED RAVE)
	SET(GFRAVE_DICTIONARY_SOURCES
	  ${CMAKE_CURRENT_SOURCE_DIR}/GFRave/include/GFRaveTrackParameters.h
	  ${CMAKE_CURRENT_SOURCE_DIR}/GFRave/include/GFRaveVertex.h
	)
	ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
	                         "${GFRAVE_DICTIONARY_SOURCES}"
	                         "${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/GFRave/include"
				 "${CMAKE_CURRENT_SOURCE_DIR}/GFRave/src/GFRaveLinkDef.h"
				 "${CMAKE_CURRENT_BINARY_DIR}/GFRaveRootDict.cc")
	SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/GFRaveRootDict.cc" "${library_sources}")
endif()


# install
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources})
INSTALL( TARGETS ${PROJECT_NAME} 
	LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
)

# link ROOT libraries
if(DEFINED RAVE)
	TARGET_LINK_LIBRARIES(
		${PROJECT_NAME}
		${ROOT_LIBS}  
  		${RAVE_LDFLAGS_STR}
  		-I${RAVE_INCLUDE_DIRS_I}
  	)
else()
	TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBS})
endif()


# the examples, build with "make tests"
ADD_CUSTOM_TARGET( tests )
MACRO( ADD_GENFIT_TEST _testname )
   ADD_EXECUTABLE( ${_testname} EXCLUDE_FROM_ALL ${ARGN} )
   ADD_DEPENDENCIES( tests  ${_testname} )
   TARGET_LINK_LIBRARIES( ${_testname} ${PROJECT_NAME}  ${ROOT_LIBS} -lGeom )
    #INSTALL( TARGETS ${_testname} DESTINATION ${EXECUTABLE_INSTALL_DIR})
ENDMACRO( ADD_GENFIT_TEST )

ADD_GENFIT_TEST( fitterTests               test/fitterTests/main.cc)
ADD_GENFIT_TEST( minimalFittingExample     test/minimalFittingExample/main.cc)
ADD_GENFIT_TEST( measurementFactoryExample test/measurementFactoryExample/main.cc)
ADD_GENFIT_TEST( streamerTest              test/streamerTest/main.cc)
ADD_GENFIT_TEST( unitTests                 test/unitTests/main.cc)
if(DEFINED RAVE)
  ADD_GENFIT_TEST( vertexingTest           test/vertexingTest/main.cc)
  ADD_GENFIT_TEST( vertexingTestRead       test/vertexingTest/read.cc)
endif()

# generate and install following configuration files
#GENERATE_PACKAGE_CONFIGURATION_FILES( genfitConfig.cmake genfitConfigVersion.cmake )
