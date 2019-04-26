######################################################################
##
## CMakeLists.txt ---
## Sample CMake project using G+Smo
##
## To ensure that G+Smo is located correctly, set the variable gismo_DIR
## to a path containing the binary and include files, eg:
##
## cmake . -Dgismo_DIR=/path/to/gismo/build
##
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0048 NEW)
endif()
if(POLICY CMP0054)# CMake 3.1
cmake_policy(SET CMP0054 NEW)
endif()
if(POLICY CMP0053)# CMake 3.1.3
cmake_policy(SET CMP0053 NEW)
endif()
if(POLICY CMP0063)# CMake 3.3
cmake_policy(SET CMP0063 NEW)
endif()

project(main)

if (NOT DEFINED CMAKE_CXX_STANDARD)
	set (CMAKE_CXX_STANDARD 14)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-enum-compare")

# G+Smo
find_package(gismo REQUIRED)
include_directories( ${GISMO_INCLUDE_DIRS})
link_directories   ( ${GISMO_LIBRARY_DIR} )
include( ${GISMO_USE_FILE})
message(STATUS "G+Smo shared library in ${GISMO_LIBRARY_DIR}")

include( ${GISMO_USE_FILE} )
set(CMAKE_MODULE_PATH ${GISMO_MODULE_PATH})
#include(gsConfig) #compiler configuration

if(NOT BASE_FOLDER)
  set(BASE_FOLDER ${CMAKE_SOURCE_DIR} CACHE STRING "BASE_FOLDER")
endif()
  add_definitions(-DBASE_FOLDER="${BASE_FOLDER}")
  message(STATUS "Paths are relative to ${BASE_FOLDER}")
 

# SoelWithDebInfi -->NDEBUG
# Remove NDEBUG flag from RelWithDebInfo builds
if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()

file(GLOB SOURCES ./gs*.cpp)
add_executable(main main.cpp ${SOURCES})
set_target_properties(main PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(main gismo)

message(STATUS "Type \"make\" to build main")
