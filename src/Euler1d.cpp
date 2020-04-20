#include "Euler1d.hpp"

Euler1dState::Euler1dState(const Mesh &mesh):
        mesh(mesh),
        u(mesh),
        pressure(mesh),
        rho(mesh),
        Energy(mesh){}

Euler1dState::Euler1dState(double left, double right, size_t points_number):
        mesh(points_number, left, right),
        u(mesh),
        pressure(mesh),
        rho(mesh),
        Energy(mesh){}

void Euler1dState::SetInitialValues(
        Euler1dState::Functor initial_u,
        Euler1dState::Functor initial_rho,
        Euler1dState::Functor initial_pressure,
        double initial_time)
{
    u.SetInitialValue(initial_u, initial_time);
    rho.SetInitialValue(initial_rho, initial_time);
    pressure.SetInitialValue(initial_pressure, initial_time);
    this->initial_time = initial_time;
}

