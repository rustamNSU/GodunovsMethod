/* *********************************************************************** */
/* this is the main function to run all tests
\* *********************************************************************** */

#include "gtest/gtest.h"

int main( int argc, char **argv ){

    /* run tests */
    ::testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS();

    return 0;

}