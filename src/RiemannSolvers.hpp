#pragma once

#include <cmath>
#include <algorithm>

struct VariableVector3 {
    double x; // mesh element center
    double q1;
    double q2;
    double q3;

    VariableVector3() = default;

    VariableVector3(double x, double q1, double q2, double q3) :
            x(x), q1(q1), q2(q2), q3(q3) {}

    friend VariableVector3 operator+(const VariableVector3 &v1, const VariableVector3 &v2);

    friend VariableVector3 operator-(const VariableVector3 &v1, const VariableVector3 &v2);

    friend VariableVector3 operator*(double scalar, const VariableVector3 &v);
};

VariableVector3 CalculateFlux(
        double velocity,
        double enthalpy,
        double gamma,
        const VariableVector3 &conservative_variable);

struct EulerVariables {
    double u;           // Velocity
    double rho;         // Density
    double pressure;
    double energy;      // Energy
    double x;           // Mesh element center
    double sound_speed; // The sound speed
    double enthalpy;    // The total specific enthalpy
    double gamma;       // The adiabatic index

    EulerVariables(
            double center_coordinate,
            double velocity,
            double density,
            double pressure,
            double gamma);

    VariableVector3 ConservativeVariable();
};


VariableVector3 RoeSolver(EulerVariables f_left, EulerVariables f_right);
