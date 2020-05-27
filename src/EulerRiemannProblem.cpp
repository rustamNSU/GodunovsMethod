#include "EulerRiemannProblem.hpp"

EulerRiemannProblem::EulerRiemannProblem() = default;

EulerRiemannProblem::EulerRiemannProblem(
        double rhol,
        double ul,
        double pl,
        double rhor,
        double ur,
        double pr,
        double gamma)
{
    setInitialValues(rhol, ul, pl, rhor, ur, pr, gamma);
}

void EulerRiemannProblem::setInitialValues(
        double rhol,
        double ul,
        double pl,
        double rhor,
        double ur,
        double pr,
        double gamma)
{
    this->rhol = rhol;
    this->ul = ul;
    this->pl = pl;
    this->rhor = rhor;
    this->ur = ur;
    this->pr = pr;
    this->gamma = gamma;

    Al = 2.0 / (gamma + 1.0) / rhol;
    Ar = 2.0 / (gamma + 1.0) / rhor;
    Bl = (gamma - 1.0) * pl / (gamma + 1.0);
    Br = (gamma - 1.0) * pr / (gamma + 1.0);
    cl = std::sqrt(gamma * pl / rhol);
    cr = std::sqrt(gamma * pr / rhor);
}

EulerVariables EulerRiemannProblem::operator()(double x, double t)
{
    if (!check_solving)
    {
        DefineRiemannWaveType();
    }

    if (x < inner_velocity * t)
    {
        /* left wave */
        if (check_wave_type_left)
        {
            /* shock wave */
            if (x < Sl * t)
            {
                return EulerVariables(x, ul, rhol, pl, gamma);
            } else{
                return EulerVariables(x, inner_velocity, inner_rhol, inner_pressure, gamma);
            }
        } else
        {
            /* rarefaction wave */
            if (x <= Shl * t)
            {
                return EulerVariables(x, ul, rhol, pl, gamma);
            } else if (x >= Stl * t)
            {
                return EulerVariables(x, inner_velocity, inner_rhol, inner_pressure, gamma);
            } else
            {
                double coef1 = 2.0 / (gamma + 1.0);
                double coef2 = (gamma - 1.0) / (gamma + 1.0);
                return EulerVariables(
                        x,
                        coef1 * (cl + 0.5 * (gamma - 1.0) * ul + x / t),
                        rhol * std::pow(coef1 + coef2 / cl * (ul - x / t), 2.0 / (gamma - 1.0)),
                        pl * std::pow(coef1 + coef2 / cl * (ul - x / t), 2.0 * gamma / (gamma - 1.0)),
                        gamma);
            }
        }
    } else
    {
        /* right wave */
        if (check_wave_type_right)
        {
            /* shock wave */
            if (x > Sr * t)
            {
                return EulerVariables(x, ur, rhor, pr, gamma);
            } else{
                return EulerVariables(x, inner_velocity, inner_rhor, inner_pressure, gamma);
            }
        } else
        {
            /* rarefaction wave */
            if (x >= Shr * t)
            {
                return EulerVariables(x, ur, rhor, pr, gamma);
            } else if (x <= Str * t)
            {
                return EulerVariables(x, inner_velocity, inner_rhor, inner_pressure, gamma);
            } else
            {
                double coef1 = 2.0 / (gamma + 1.0);
                double coef2 = (gamma - 1.0) / (gamma + 1.0);
                return EulerVariables(
                        x,
                        coef1 * (-cr + 0.5 * (gamma - 1.0) * ur + x / t),
                        rhor * std::pow(coef1 - coef2 / cr * (ur - x / t), 2.0 / (gamma - 1.0)),
                        pr * std::pow(coef1 + coef2 / cr * (ur - x / t), 2.0 * gamma / (gamma - 1.0)),
                        gamma);
            }
        }
    }
}

double EulerRiemannProblem::FunctionLeft(double p)
{
    if (p > pl) // shock wave
    {
        return (p - pl) * std::sqrt(Al / (p + Bl));
    } else // rarefaction wave
    {
        return 2.0 * cl / (gamma - 1.0) * (std::pow(p / pl, 0.5 * (gamma - 1.0) / gamma) - 1.0);
    }
}

double EulerRiemannProblem::FunctionRight(double p)
{
    if (p > pr) // shock wave
    {
        return (p - pr) * std::sqrt(Ar / (p + Br));
    } else // rarefaction wave
    {
        return 2.0 * cr / (gamma - 1.0) * (std::pow(p / pr, 0.5 * (gamma - 1.0) / gamma) - 1.0);
    }
}

double EulerRiemannProblem::DerivativeFunctionLeft(double p)
{
    if (p > pl) // shock wave
    {
        return (1.0 - 0.5 * (p - pl) / (Bl + p)) * std::sqrt(Al / (p + Bl));
    } else // rarefaction wave
    {
        return std::pow(p / pl, -0.5 * (gamma + 1.0) / gamma) / rhol / cl;
    }
}

double EulerRiemannProblem::DerivativeFunctionRight(double p)
{
    if (p > pr) // shock wave
    {
        return (1.0 - 0.5 * (p - pr) / (Br + p)) * std::sqrt(Ar / (p + Br));
    } else // rarefaction wave
    {
        return std::pow(p / pr, -0.5 * (gamma + 1.0) / gamma) / rhor / cr;
    }
}

double EulerRiemannProblem::p0TR()
{
    double numerator = cl + cr - 0.5 * (gamma - 1.0) * (ur - ul);
    double denominator = std::pow(cl / pl, 0.5 * (gamma - 1.0) / gamma) +
                         std::pow(cr / pr, 0.5 * (gamma - 1.0) / gamma);
    return std::pow(numerator / denominator, 2 * gamma / (gamma - 1.0));
}

void EulerRiemannProblem::FindInnerPressure(double p0, double tolerance)
{
    while (true)
    {
        double p1 = p0 - (FunctionLeft(p0) + FunctionRight(p0) + ur - ul) /
                         (DerivativeFunctionLeft(p0) + DerivativeFunctionRight(p0));
        if (0.5 * std::abs((p1 - p0) / (p1 + p0)) < tolerance)
        {
            this->inner_pressure = p1;
            return;
        }
        p0 = p1;
    }
}

void EulerRiemannProblem::CalculateInnerVelocity()
{
    this->inner_velocity = 0.5 * (ul + ur + FunctionRight(inner_pressure) - FunctionLeft(inner_pressure));
}

void EulerRiemannProblem::DefineRiemannWaveType()
{
    FindInnerPressure(p0TR());
    CalculateInnerVelocity();
    check_solving = true; // Calculate all parameters only once

    if (inner_pressure > pl)
    {
        /* shock wave */
        this->check_wave_type_left = true;
        double coef1 = inner_pressure / pl;
        double coef2 = (gamma - 1.0) / (gamma + 1.0);

        this->Sl = ul - cl * std::sqrt(0.5 * (gamma + 1.0) / gamma * (coef1 + coef2));
        this->inner_rhol = rhol * (coef1 + coef2) / (coef1 * coef2 + 1.0);
    } else
    {
        /* rarefaction wave */
        this->check_wave_type_left = false;
        double coef1 = inner_pressure / pl;
        double coef2 = (gamma - 1.0) / (gamma + 1.0);

        this->inner_rhol = rhol * std::pow(coef1, 1.0 / gamma);
        this->Shl = ul - cl;
        this->Stl = inner_velocity - cl * std::pow(coef1, 0.5 * (gamma - 1.0) / gamma);
    }

    if (inner_pressure > pr)
    {
        /* shock wave */
        this->check_wave_type_right = true;
        double coef1 = inner_pressure / pr;
        double coef2 = (gamma - 1.0) / (gamma + 1.0);

        this->Sr = ur + cr * std::sqrt(0.5 * (gamma + 1.0) / gamma * (coef1 + coef2));
        this->inner_rhor = rhor * (coef1 + coef2) / (coef1 * coef2 + 1.0);
    } else
    {
        /* rarefaction wave */
        this->check_wave_type_right = false;
        double coef1 = inner_pressure / pr;
        double coef2 = (gamma - 1.0) / (gamma + 1.0);

        this->inner_rhor = rhor * std::pow(coef1, 1.0 / gamma);
        this->Shr = ur + cr;
        this->Str = inner_velocity + cr * std::pow(coef1, 0.5 * (gamma - 1.0) / gamma);
    }
}