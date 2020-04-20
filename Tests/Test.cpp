// #define _USE_MATH_DEFINES
#include "gtest/gtest.h"
#include "Test.hpp"
#include "Euler1d.hpp"

TEST(Euler1dTest, InitialTest)
{
    Euler1dState model(0.0, 1.0, 40);
    model.SetInitialValues(
            [](double x){return x;},
            [](double x){return x;},
            [](double x){return x;}
            );
}
