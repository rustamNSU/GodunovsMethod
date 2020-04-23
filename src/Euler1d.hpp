#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "Mesh.hpp"
#include "RiemannSolvers.hpp"
#include "EulerRiemannProblem.hpp"

class Euler1dState {
private:
    Function u; // The flow velocity
    Function rho; // The density
    Function pressure; // The pressure
    Function energy; // The Energy
    Function sound_speed; // The sound speed
    Mesh mesh;

    double gamma = 1.4; // The adiabatic index
    double initial_time = 0.0;
    double time_step;

    /* Addition */
    using Functor = double (*)(double);
    using RiemannSolver = VariableVector3(*)(const EulerVariables &, const EulerVariables &);

public:
    Euler1dState() = delete;

    Euler1dState(const Mesh &mesh);

    Euler1dState(double left, double right, size_t points_number);

    void SetInitialValues(
            Functor initial_u,
            Functor initial_rho,
            Functor initial_pressure,
            double initial_time = 0.0);

    Function GetVelocity() const;

    Function GetDensity() const;

    Function GetPressure() const;

    void SetTimeStep(double CFL = 0.5);

    void CalculateNextLayer(RiemannSolver method);

    void CalculateNextLayer(RiemannSolver method, double time_step);

    void Calculate(RiemannSolver method, double end_time);

    std::vector<EulerVariables> GenerateEulerVariable();

    void AddNextLayer(std::vector<EulerVariables> euler_variable, double time_step);
};