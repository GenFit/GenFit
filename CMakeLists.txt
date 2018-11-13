CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
if(POLICY CMP0042)
   # policies not known to CMake 2.8 ...
   cmake_policy(SET CMP0042 NEW) # Mac OS -rpath behavior
endif(POLICY CMP0042)

# declare the project name
PROJECT(genfit2)

# set project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 2 )
SET( ${PROJECT_NAME}_VERSION_MINOR 2 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

include(GNUInstallDirs)

# install destinations can be passed via the command line:
# cmake -DLIBRARY_OUTPUT_DIRECTORY:PATH=.... -DINCLUDE_OUTPUT_DIRECTORY:PATH=....
# otherwise, set default install destinations
IF(NOT DEFINED LIBRARY_OUTPUT_DIRECTORY)
	set( LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_LIBDIR} )
ENDIF()
IF(NOT DEFINED INCLUDE_OUTPUT_DIRECTORY)
	set( INCLUDE_OUTPUT_DIRECTORY ${CMAKE_INSTALL_INCLUDEDIR} )
ENDIF()

# define path to macros and settings
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
# include the macros and settings
INCLUDE(settings)
INCLUDE(macros)
INCLUDE(genfit)

FIND_PACKAGE(ROOT REQUIRED COMPONENTS Geom)
IF(${ROOT_found_version} LESS 60000)
	MESSAGE(FATAL_ERROR "GenFit requires a ROOT version > 6.00.00")
ENDIF()

FIND_PACKAGE(Eigen3 REQUIRED)

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

INSTALL(DIRECTORY ${GF_INC_DIRS}
        DESTINATION ${INCLUDE_OUTPUT_DIRECTORY}
        PATTERN ".svn" EXCLUDE)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CPP_STANDARD c++11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=${CPP_STANDARD} -pthread -Wpedantic -Wall -Wextra \
 -Wno-long-long -Wshadow -Werror=overloaded-virtual -Wunused-parameter -Waddress")

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O3 -g")


INCLUDE_DIRECTORIES(
		BEFORE
		${ROOT_INCLUDE_DIRS}
		SYSTEM ${EIGEN3_INCLUDE_DIR}
		${GF_INC_DIRS}
)

FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
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
SET(CORE_DICTIONARY_SOURCES
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsFinitePlane.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsFitterInfo.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsHMatrix.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/AbsTrackRep.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/DetPlane.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/FitStatus.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/IO.h
		${CMAKE_CURRENT_SOURCE_DIR}/core/include/Material.h
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
ROOT_GENERATE_DICTIONARY(
		"${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
		"${CORE_DICTIONARY_SOURCES}"
		"${CMAKE_CURRENT_SOURCE_DIR}/core/include"
		"${CMAKE_CURRENT_SOURCE_DIR}/core/src/coreLinkDef.h"
		"${CMAKE_CURRENT_BINARY_DIR}/coreRootDict.cc"
)
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/coreRootDict.cc" "${library_sources}")

SET(EVENTDISPLAY_DICTIONARY_SOURCES
		${CMAKE_CURRENT_SOURCE_DIR}/eventDisplay/include/EventDisplay.h
)
ROOT_GENERATE_DICTIONARY(
		"${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
		"${EVENTDISPLAY_DICTIONARY_SOURCES}"
		"${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/fitters/include;${CMAKE_CURRENT_SOURCE_DIR}/eventDisplay/include"
		"${CMAKE_CURRENT_SOURCE_DIR}/eventDisplay/src/eventDisplayLinkDef.h"
		"${CMAKE_CURRENT_BINARY_DIR}/eventDisplayRootDict.cc"
)
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
ROOT_GENERATE_DICTIONARY(
		"${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
		"${FITTERS_DICTIONARY_SOURCES}"
		"${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/fitters/include"
		"${CMAKE_CURRENT_SOURCE_DIR}/fitters/src/fittersLinkDef.h"
		"${CMAKE_CURRENT_BINARY_DIR}/fittersRootDict.cc"
)
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/fittersRootDict.cc" "${library_sources}")

SET(GBL_DICTIONARY_SOURCES
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GFGbl.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblFitter.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/ICalibrationParametersDerivatives.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblFitterInfo.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblFitStatus.h
		${CMAKE_CURRENT_SOURCE_DIR}/GBL/include/GblData.h
)
ROOT_GENERATE_DICTIONARY(
		"${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
        "${GBL_DICTIONARY_SOURCES}"
		"${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/GBL/include"
		"${CMAKE_CURRENT_SOURCE_DIR}/GBL/src/GBLLinkDef.h"
		"${CMAKE_CURRENT_BINARY_DIR}/GBLRootDict.cc"
)
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
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WireMeasurementNew.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WirePointMeasurement.h
		${CMAKE_CURRENT_SOURCE_DIR}/measurements/include/WireTrackCandHit.h
)
ROOT_GENERATE_DICTIONARY(
		"${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
		"${MEASUREMENTS_DICTIONARY_SOURCES}"
        "${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/measurements/include"
		"${CMAKE_CURRENT_SOURCE_DIR}/measurements/src/measurementsLinkDef.h"
		"${CMAKE_CURRENT_BINARY_DIR}/measurementsRootDict.cc"
)
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/measurementsRootDict.cc" "${library_sources}")

SET(TRACKREPS_DICTIONARY_SOURCES
		${CMAKE_CURRENT_SOURCE_DIR}/trackReps/include/RKTrackRep.h
)
ROOT_GENERATE_DICTIONARY(
		"${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
		"${TRACKREPS_DICTIONARY_SOURCES}"
		"${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/trackReps/include"
		"${CMAKE_CURRENT_SOURCE_DIR}/trackReps/src/trackRepsLinkDef.h"
		"${CMAKE_CURRENT_BINARY_DIR}/trackRepsRootDict.cc"
)
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/trackRepsRootDict.cc" "${library_sources}")

SET(UTILITIES_DICTIONARY_SOURCES
		${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/HelixTrackModel.h
		${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/MeasurementCreator.h
		${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/mySpacepointDetectorHit.h
		${CMAKE_CURRENT_SOURCE_DIR}/utilities/include/mySpacepointMeasurement.h
)
ROOT_GENERATE_DICTIONARY(
		"${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
		"${UTILITIES_DICTIONARY_SOURCES}"
		"${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/utilities/include;${CMAKE_CURRENT_SOURCE_DIR}/measurements/include"
		"${CMAKE_CURRENT_SOURCE_DIR}/utilities/src/utilitiesLinkDef.h"
		"${CMAKE_CURRENT_BINARY_DIR}/utilitiesRootDict.cc"
)
SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/utilitiesRootDict.cc" "${library_sources}")

if(DEFINED RAVE)
	SET(GFRAVE_DICTIONARY_SOURCES
			${CMAKE_CURRENT_SOURCE_DIR}/GFRave/include/GFRaveTrackParameters.h
			${CMAKE_CURRENT_SOURCE_DIR}/GFRave/include/GFRaveVertex.h
	)
	ROOT_GENERATE_DICTIONARY(
			"${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
			"${GFRAVE_DICTIONARY_SOURCES}"
			"${CMAKE_CURRENT_SOURCE_DIR}/core/include;${CMAKE_CURRENT_SOURCE_DIR}/GFRave/include"
			"${CMAKE_CURRENT_SOURCE_DIR}/GFRave/src/GFRaveLinkDef.h"
			"${CMAKE_CURRENT_BINARY_DIR}/GFRaveRootDict.cc"
	)
	SET(library_sources "${CMAKE_CURRENT_BINARY_DIR}/GFRaveRootDict.cc" "${library_sources}")
endif()


# install
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources})
INSTALL(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

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


ADD_CUSTOM_TARGET( tests )

ADD_GENFIT_TEST( fitterTests               test/fitterTests/main.cc)
ADD_GENFIT_TEST( minimalFittingExample     test/minimalFittingExample/main.cc)
ADD_GENFIT_TEST( minimalFittingExample2    test/minimalFittingExample/main2.cc)
ADD_GENFIT_TEST( measurementFactoryExample test/measurementFactoryExample/main.cc)
ADD_GENFIT_TEST( streamerTest              test/streamerTest/main.cc)
ADD_GENFIT_TEST( unitTests                 test/unitTests/main.cc)
IF(DEFINED RAVE)
  ADD_GENFIT_TEST( vertexingTest           test/vertexingTest/main.cc)
  ADD_GENFIT_TEST( vertexingTestRead       test/vertexingTest/read.cc)
ENDIF()


FIND_PACKAGE(GTest)
IF(GTEST_FOUND)
	ADD_EXECUTABLE(gtests
			gtest/TestConstField.cpp
			gtest/TestMaterialEffects.cpp
			gtest/TestRKTrackRep.cpp
			gtest/TestRootEigenTransformations.cpp
			gtest/TestRKMatrixEigenTransformations.cpp
			gtest/TestUnits.cpp
			gtest/TestMaterial.cpp
			)
	TARGET_LINK_LIBRARIES(gtests ${GTEST_BOTH_LIBRARIES} ${ROOT_LIBS} ${PROJECT_NAME})  # gtest gtest_main
	MESSAGE(STATUS  ${GTEST_INCLUDE_DIRS})
	TARGET_INCLUDE_DIRECTORIES(gtests PUBLIC ${GTEST_INCLUDE_DIRS})
	INSTALL(TARGETS gtests DESTINATION bin)
ELSE()
	MESSAGE(STATUS "GTest not installed or found -- skip building unittests.")
ENDIF()

# generate and install following configuration files
#GENERATE_PACKAGE_CONFIGURATION_FILES( genfitConfig.cmake genfitConfigVersion.cmake )
