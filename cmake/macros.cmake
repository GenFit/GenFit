# collection of helpful macros for cmake
# created by J. Engels, DESY and Ch. Rosemann, DESY

EXECUTE_PROCESS( COMMAND "svnversion"
        OUTPUT_VARIABLE GLOBAL_SVN_REVISION
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        RESULT_VARIABLE _exit_code
        )
  IF( NOT _exit_code EQUAL 0 )
        MESSAGE( STATUS "Couldn't retrieve a version number from SVN to set in documentation.")
        SET(GLOBAL_SVN_REVISION "NoNumberAvailable")
    ENDIF()

MESSAGE(STATUS "Found global svn revision to be ${GLOBAL_SVN_REVISION}.")


# create symbolic lib target for calling library targets
MACRO( ADD_SHARED_LIBRARY _name )
    ADD_LIBRARY( ${_name} SHARED ${ARGN} )
    
    # change lib_target properties
    SET_TARGET_PROPERTIES( ${_name} PROPERTIES
        # create *nix style library versions + symbolic links
        VERSION ${${PROJECT_NAME}_VERSION}
        SOVERSION ${${PROJECT_NAME}_SOVERSION}
    )
ENDMACRO( ADD_SHARED_LIBRARY )


# in order to include cmake projects into other projects config files must exist
# helper macro for generating project configuration file
MACRO( GENERATE_PACKAGE_CONFIGURATION_FILES )

    FOREACH( arg ${ARGN} )
        IF( ${arg} MATCHES "Config.cmake" )
            IF( EXISTS "${PROJECT_SOURCE_DIR}/cmake/${arg}.in" )
                CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/${arg}.in"
                                "${PROJECT_BINARY_DIR}/${arg}" @ONLY
                )
                INSTALL( FILES "${PROJECT_BINARY_DIR}/${arg}" DESTINATION . )
             ENDIF()
        ENDIF()


        IF( ${arg} MATCHES "ConfigVersion.cmake" )
            # version configuration file
            IF( EXISTS "${PROJECT_SOURCE_DIR}/cmake/${arg}.in" )
                CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/${arg}.in"
                                "${PROJECT_BINARY_DIR}/${arg}" @ONLY
                )
                INSTALL( FILES "${PROJECT_BINARY_DIR}/${arg}" DESTINATION . )
            ENDIF( EXISTS "${PROJECT_SOURCE_DIR}/cmake/${arg}.in" )
        ENDIF()

        IF( ${arg} MATCHES "LibDeps.cmake" )
            EXPORT_LIBRARY_DEPENDENCIES( "${arg}" )
            INSTALL( FILES "${PROJECT_BINARY_DIR}/${arg}" DESTINATION lib/cmake )
        ENDIF()

    ENDFOREACH()

ENDMACRO( GENERATE_PACKAGE_CONFIGURATION_FILES )
