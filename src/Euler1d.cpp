#include "Euler1d.hpp"
#include <cmath>
#include <algorithm>

Euler1dState::Euler1dState(const Mesh &mesh) :
        mesh(mesh),
        u(mesh),
        pressure(mesh),
        rho(mesh),
        energy(mesh),
        sound_speed(mesh)
{
}

Euler1dState::Euler1dState(double left, double right, size_t points_number) :
        mesh(points_number, left, right),
        u(mesh),
        pressure(mesh),
        rho(mesh),
        energy(mesh),
        sound_speed(mesh)
{
}

void Euler1dState::SetInitialValues(
        Euler1dState::Functor initial_u,
        Euler1dState::Functor initial_rho,
        Euler1dState::Functor initial_pressure,
        double initial_time)
{
    u.SetInitialValue(initial_u, initial_time);
    rho.SetInitialValue(initial_rho, initial_time);
    pressure.SetInitialValue(initial_pressure, initial_time);

    auto slice_function_u = u.GetLastLayer();
    auto slice_function_rho = rho.GetLastLayer();
    auto slice_function_pressure = pressure.GetLastLayer();
    auto initial_energy = (1.0 / (gamma - 1.0)) * slice_function_pressure +
                          0.5 * slice_function_rho * slice_function_u * slice_function_u;

    auto initial_sound_speed = gamma * slice_function_pressure / slice_function_rho;
    initial_sound_speed.ApplyForeach(std::sqrt);
    energy.SetInitialValue(initial_energy, initial_time);
    sound_speed.SetInitialValue(initial_sound_speed, initial_time);
    this->initial_time = initial_time;
    SetTimeStep();
}

Function Euler1dState::GetVelocity() const
{
    return this->u;
}

Function Euler1dState::GetDensity() const
{
    return this->rho;
}

Function Euler1dState::GetPressure() const
{
    return this->pressure;
}

void Euler1dState::SetTimeStep(double CFL)
{
    auto eigen_value1 = u.GetLastLayer() - sound_speed.GetLastLayer();
    auto eigen_value3 = u.GetLastLayer() + sound_speed.GetLastLayer();
    eigen_value1.ApplyForeach(std::abs);
    eigen_value3.ApplyForeach(std::abs);

    double max_wave_speed = 0.0;
    for (size_t i = 0; i < mesh.GetSize(); ++i)
    {
        double wave_speed = std::max(eigen_value1[i], eigen_value3[i]);
        if (wave_speed > max_wave_speed) max_wave_speed = wave_speed;
    }
    time_step = CFL * mesh.GetStep() / max_wave_speed;
}

void Euler1dState::CalculateNextLayer(Euler1dState::RiemannSolver method)
{
    this->CalculateNextLayer(method, this->time_step);
}

void Euler1dState::CalculateNextLayer(Euler1dState::RiemannSolver method, double time_step)
{
    std::vector<VariableVector3> inner_flux;
    inner_flux.reserve(mesh.GetSize() - 1);
    auto EulerVariable = this->GenerateEulerVariable();
    for (size_t i = 0; i < mesh.GetSize() - 1; ++i)
    {
        inner_flux.push_back(
                method(EulerVariable[i], EulerVariable[i+1]));
    }
    /*
     * We consider that the wave doesn't reach the domain boundary,
     * therefore we transfer the value in the boundary cells to the next layer
     */
    std::vector<VariableVector3> ConservativeVariable;
    ConservativeVariable.reserve(EulerVariable.size());
    for (size_t i = 0; i < EulerVariable.size(); ++i)
    {
        ConservativeVariable.push_back(EulerVariable[i].ConservativeVariable());
    }

    /* The integral form approximation (Godunov step) */
    for (size_t i = 0; i < inner_flux.size() - 1; ++i)
    {
        ConservativeVariable[i+1] = ConservativeVariable[i+1] - (time_step / mesh.GetStep()) *
                                                                (inner_flux[i+1] - inner_flux[i]);
    }

    /* Recovering equation parameters from conservative variables */
    for (size_t i = 0; i < EulerVariable.size(); ++i)
    {
        EulerVariable[i].RestoreFromConservativeVariable(ConservativeVariable[i]);
    }

    AddNextLayer(EulerVariable, time_step);
}

void Euler1dState::Calculate(RiemannSolver method, double end_time)
{
    double time = u.GetLastTimeStep() + time_step;
    while (time < end_time)
    {
        CalculateNextLayer(method);
        time += time_step;
    }
}

std::vector<EulerVariables> Euler1dState::GenerateEulerVariable()
{
    std::vector<EulerVariables> result;
    result.reserve(mesh.GetSize());
    for (size_t i = 0; i < mesh.GetSize(); ++i)
    {
        result.emplace_back(
                mesh[i],
                u.GetLastLayer()[i],
                rho.GetLastLayer()[i],
                pressure.GetLastLayer()[i],
                gamma
                );
    }
    return result;
}

void Euler1dState::AddNextLayer(std::vector<EulerVariables> euler_variable, double time_step)
{
    std::vector<double> u_next;
    std::vector<double> rho_next;
    std::vector<double> pressure_next;
    std::vector<double> energy_next;
    std::vector<double> sound_speed_next;

    u_next.reserve(euler_variable.size());
    rho_next.reserve(euler_variable.size());
    pressure_next.reserve(euler_variable.size());
    energy_next.reserve(euler_variable.size());
    sound_speed_next.reserve(euler_variable.size());

    for (auto &elem : euler_variable)
    {
        u_next.push_back(elem.u);
        rho_next.push_back(elem.rho);
        pressure_next.push_back(elem.pressure);
        energy_next.push_back(elem.energy);
        sound_speed_next.push_back(elem.sound_speed);
    }

    u.AddLayer(time_step, SliceFunction(std::move(u_next)));
    rho.AddLayer(time_step, SliceFunction(std::move(rho_next)));
    pressure.AddLayer(time_step, SliceFunction(std::move(pressure_next)));
    energy.AddLayer(time_step, SliceFunction(std::move(energy_next)));
    sound_speed.AddLayer(time_step, SliceFunction(std::move(sound_speed_next)));
}

