#pragma once

#include <vector>
#include <cmath>
#include "Mesh.hpp"
#include "RiemannSolvers.hpp"


class Euler1dState {
private:
    Function u; // The flow velocity
    Function rho; // The density
    Function pressure; // The pressure
    Function Energy; // The Energy
    Mesh mesh;

    double gamma = 1.4; // The adiabatic index
    double initial_time = 0.0;
    double time_step;

    /* Addition function */
    using Functor = double (*)(double);

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

    template<typename Method>
    void CalculateNextLayer(Method* method);

    template<typename Method>
    void CalculateNextLayer(Method* method, double time_step);
};