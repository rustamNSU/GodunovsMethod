cmake_minimum_required (VERSION 3.8)
project (GodunovMethod)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Include cmake module ExternalProject needed to build external dependencies
include(ExternalProject)    

# Enable 11 C++ standard
set(CMAKE_CXX_STANDARD 17) 
set(CMAKE_C_FLAGS -m64)
set(CMAKE_CXX_FLAGS -m64)

if(NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif()
if(NOT CMAKE_CONFIGURATION_TYPES)
    set (CMAKE_CONFIGURATION_TYPES Release)
endif()

# add_compile_options(-Wall -Wextra -pedantic -Werror)

# INSTALL_PREFIX MATCHES "^/usr/local" OR CMAKE_INSTALL_PREFIX MATCHES "^C:/Program Files" )
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/Installation)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


# RPath folder
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

include(${CMAKE_CURRENT_SOURCE_DIR}/externalGoogleTest.cmake)
add_subdirectory (src)
add_subdirectory (Tests) 