cmake_minimum_required(VERSION 3.15)
project(ProteinElectricField)

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release ... FORCE )
endif()

set(CMAKE_CXX_STANDARD 20)

# Set the binary output
if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/release)
endif()

#-----------------------------------------------[INCLUDE GIT  SUBMODULES]-----------------------------------------------

add_subdirectory(libs/spdlog)
add_subdirectory(libs/cxxopts)

#----------------------------------------------[INCLUDE COMPILER WARNINGS]----------------------------------------------

include(cmake/compilerWarnings.cmake)
include(cmake/sanitizers.cmake)
include(cmake/StandardCmake.cmake)

add_subdirectory(src)
