#pragma once

#include "RiemannSolvers.hpp"
#include "Mesh.hpp"
#include <cmath>

class EulerRiemannProblem {
    double ul; // The initial velocity x < 0
    double ur; // The initial velocity x > 0
    double rhol; // The initial density x < 0
    double rhor; // The initial density x > 0
    double pl; // The initial pressure x < 0
    double pr; // The initial pressure x > 0

    double gamma; // The adiabatic coefficient

    /* Addition coefficient */
    double Al;
    double Ar;
    double Bl;
    double Br;
    double cl; // The sound speed (left)
    double cr; // The sound speed (right)

    /* Star-region solution (inner) */
    bool check_solving = false;
    double inner_pressure; // Find from solving equation
    double inner_velocity; // Define after finding inner_pressure

    bool check_wave_type_left; // true == shock, false == rarefaction
    double inner_rhol;
    double Sl; // The left shock wave speed
    double Shl; // The left rarefaction Head-wave speed
    double Stl; // The left rarefaction Tail-wave speed

    bool check_wave_type_right; // true == shock, false == rarefaction
    double inner_rhor;
    double Sr; // The right shock wave speed
    double Shr; // The right rarefaction Head-wave speed
    double Str; // The right rarefaction Tail-wave speed

public:
    EulerRiemannProblem(
            double rhol, // The initial density x < 0
            double ul, // The initial velocity x < 0
            double pl, // The initial pressure x < 0
            double rhor, // The initial density x > 0
            double ur, // The initial velocity x > 0
            double pr, // The initial pressure x > 0
            double gamma = 1.4
            );

    EulerVariables operator()(double x, double t);

    double FunctionLeft(double p);

    double FunctionRight(double p);

    double DerivativeFunctionLeft(double p);

    double DerivativeFunctionRight(double p);

    double p0TR();

    void FindInnerPressure(double p0, double tolerance = 1e-6);

    void CalculateInnerVelocity();

    void DefineRiemannWaveType();
};
