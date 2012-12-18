cmake_minimum_required(VERSION 2.8)

# The SIMPATH variable is defined in the PANDAroot environment.
# It is used here to determine which CMakeLists to use.
IF(DEFINED ENV{SIMPATH})
   include(CMakeLists_panda.txt)
ELSE(DEFINED ENV{SIMPATH})
   include(CMakeLists_standalone.txt)
ENDIF(DEFINED ENV{SIMPATH})
