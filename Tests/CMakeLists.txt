cmake_minimum_required (VERSION 3.8)

project (GodunovMethodTest)

file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})  

# Compile executable from source files in SOURCES variable
add_executable(GodunovMethodTest ${SOURCES})         

target_link_libraries(GodunovMethodTest GodunovMethod googletest)

install (TARGETS GodunovMethodTest DESTINATION bin)
