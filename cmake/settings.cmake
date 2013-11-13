# general settings for a cmake project
# less clutter in the real CMakeLists.txt file


# library *nix style versioning
SET( ${PROJECT_NAME}_SOVERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}" )
SET( ${PROJECT_NAME}_VERSION   "${${PROJECT_NAME}_SOVERSION}.${${PROJECT_NAME}_VERSION_PATCH}" )



# define output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH )
MARK_AS_ADVANCED( LIBRARY_OUTPUT_PATH )

# what happens when `make install` is invoked:
# set default install prefix to project root directory
# instead of the cmake default /usr/local
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" )
    SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()

# write this variable to cache
SET( CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Where to install ${PROJECT_NAME}" FORCE )

