// #define _USE_MATH_DEFINES
#include "gtest/gtest.h"
#include "Test.hpp"
#include "Euler1d.hpp"

#include <iostream>

TEST(Euler1dTest, InitialTest)
{
    Euler1dState model(0.0, 1.0, 10);
    model.SetInitialValues(
            [](double x){return 1.0;},
            [](double x){return 1.0;},
            [](double x){return 1.0;}
            );
    model.CalculateNextLayer(RoeSolver);
    Function u = model.GetVelocity();
    Function rho = model.GetDensity();
    Function pressure = model.GetPressure();
}
