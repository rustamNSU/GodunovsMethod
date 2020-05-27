// #define _USE_MATH_DEFINES
#include "gtest/gtest.h"
#include "Test.hpp"
#include "Euler1d.hpp"
#include "EulerRiemannProblem.hpp"

#include <iostream>

using namespace std;

TEST(Euler1dTest, InitialTest)
{
    Euler1dState model(0.0, 1.0, 10);
    model.SetInitialValues(
            [](double x) { return 1.0; },
            [](double x) { return 1.0; },
            [](double x) { return 1.0; }
    );
    model.CalculateNextLayer(RoeSolver);
    Function u = model.GetVelocity();
    Function rho = model.GetDensity();
    Function pressure = model.GetPressure();
}

TEST(Euler1dTest, SodShockTube)
{
    auto SodInitialVelocity = [](double x) { return 0.0; };
    auto SodInitialDensity = [](double x)
    {
        if (x < 0.5) {return 1.0; }
        else {return 0.125; }
    };
    auto SodInitialPressure = [](double x)
    {
        if (x < 0.5) {return 10.0; }
        else {return 0.1; }
    };

    EulerRiemannProblem exact_solution(1.0, 0.0, 10.0,0.125, 0.0, 0.1);
    exact_solution.DefineRiemannWaveType();
    double end_time = 0.2;
    size_t N = 100;

    Euler1dState model(0.0, 1.0, N);
    model.SetInitialValues(
            SodInitialVelocity,
            SodInitialDensity,
            SodInitialPressure
            );
    model.Calculate(RoeSolver, end_time);

    auto u = model.GetVelocity();
    auto rho = model.GetDensity();
    auto pressure = model.GetPressure();

    cout << "Number of layers = " << model.GetLayersNumber() << endl;
    cout << static_cast<int>(model.GetLayersNumber() / 2) << " - layer time = "
         << model.GetLayerTime(model.GetLayersNumber() / 2) << endl;

    cout << "exact Sod:" << endl;
    cout << "  exact_solution(0.1, 1.0) = " << exact_solution(0.1, 1.0).rho << endl;
}
