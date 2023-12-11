MACRO( ADD_GENFIT_TEST _testname )
    # If we have tests enabled build all tests directly, otherwise build them
    # only via the tests target
    IF(BUILD_TESTING)
      ADD_EXECUTABLE( ${_testname} ${ARGN} )
    ELSE()
      ADD_EXECUTABLE( ${_testname} EXCLUDE_FROM_ALL ${ARGN} )
    ENDIF()
    ADD_DEPENDENCIES( tests  ${_testname} )
    TARGET_LINK_LIBRARIES( ${_testname} ${PROJECT_NAME}  ROOT::Core ROOT::Geom ROOT::Physics)
    #INSTALL( TARGETS ${_testname} DESTINATION ${EXECUTABLE_INSTALL_DIR})
ENDMACRO( ADD_GENFIT_TEST )
