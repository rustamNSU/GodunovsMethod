// #define _USE_MATH_DEFINES
#include "gtest/gtest.h"
#include "Test.hpp"
#include "Euler1d.hpp"

#include <iostream>

TEST(Euler1dTest, InitialTest)
{
    Euler1dState model(0.0, 1.0, 40);
    model.SetInitialValues(
            [](double x){return x;},
            [](double x){return x;},
            [](double x){return x;}
            );
    Function u = model.GetVelocity();
    Function rho = model.GetDensity();
    Function pressure = model.GetPressure();
    std::cout << "sizeof(Flux3) = " << sizeof(VariableVector3)
    << ", 3 * sizeof(double) = " << 4 * sizeof(double) << std::endl;
}
