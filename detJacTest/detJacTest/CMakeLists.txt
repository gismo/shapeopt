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

set (CMAKE_CXX_STANDARD 11)

# G+Smo
find_package(gismo REQUIRED)
include_directories( ${GISMO_INCLUDE_DIRS})
link_directories   ( ${GISMO_LIBRARY_DIR} )
include( ${GISMO_USE_FILE})
message(STATUS "G+Smo shared library in ${GISMO_LIBRARY_DIR}")

include( ${GISMO_USE_FILE} )
set(CMAKE_MODULE_PATH ${GISMO_MODULE_PATH})
#include(gsConfig) #compiler configuration

# Source file
file(GLOB SOURCES ./*.cpp)
add_executable(main ${SOURCES})
set_target_properties(main PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(main gismo)

message(STATUS "Type \"make\" to build main")
