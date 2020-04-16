set(GOOGLETEST_INSTALL_DIR ${CMAKE_BINARY_DIR}/googletestInstallation)

# Command to compile googletest submodule source
# SOURCE_DIR - Path to googletest source
# INSTALL_DIR - Path where to install googletest library (This variable will be used in project subdirectories)
# BINARY_DIR - Path to folder for out-of-source build of googletest library
# CMAKE_ARGS - Actual designation of install folder
ExternalProject_Add(googletest_proj
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/googletest
INSTALL_DIR ${GOOGLETEST_INSTALL_DIR} 
BINARY_DIR ${CMAKE_BINARY_DIR}/googletestBuild
BUILD_ALWAYS 1
CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DINSTALL_GTEST=1 -Dgtest_force_shared_crt=ON -DCMAKE_INSTALL_BINDIR=${GOOGLETEST_INSTALL_DIR} -DCMAKE_INSTALL_LIBDIR=${GOOGLETEST_INSTALL_DIR} -DCMAKE_INSTALL_INCLUDEDIR=${GOOGLETEST_INSTALL_DIR} -DCMAKE_INSTALL_PREFIX=${GOOGLETEST_INSTALL_DIR} -DCMAKE_INSTALL_MESSAGE=LAZY
BUILD_BYPRODUCTS ${GOOGLETEST_INSTALL_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX} ${GOOGLETEST_INSTALL_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}
)

# dummy library target
# 
add_library(googletest INTERFACE)
add_dependencies(googletest googletest_proj)

target_include_directories(googletest INTERFACE ${GOOGLETEST_INSTALL_DIR})

target_link_libraries(googletest INTERFACE ${GOOGLETEST_INSTALL_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX} 
                                           ${GOOGLETEST_INSTALL_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX} 
                                           pthread )
