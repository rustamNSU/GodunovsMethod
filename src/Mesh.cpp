#include "Mesh.hpp"

Mesh::Mesh() = default;

Mesh::Mesh(size_t points_number, double left, double right)
{
    this->SetMesh(points_number, left, right);
}

void Mesh::SetMesh(size_t points_number, double left, double right)
{
    this->mesh.resize(points_number);

    double step = (right - left) / points_number;
    this->step_x = step;
    double element_center = left + step / 2.0;
    for (auto &elem : this->mesh)
    {
        elem = element_center;
        element_center += step;
    }
}

size_t Mesh::GetSize() const
{
    return mesh.size();
}

double Mesh::GetStep() const
{
    return step_x;
}

double Mesh::operator[](size_t index) const
{
    return mesh[index];
}

double &Mesh::operator[](size_t index)
{
    return mesh[index];
}

SliceFunction::SliceFunction() = default;

SliceFunction::SliceFunction(std::vector<double> &&value) :
    value(std::move(value)) {}

SliceFunction::SliceFunction(const Mesh &mesh, double(* f)(double))
{
    this->SetValue(mesh, f);
}

void SliceFunction::SetValue(const Mesh &mesh, double(* f)(double))
{
    int index = 0;
    std::vector<double> result(mesh.GetSize());
    for (auto &elem : result)
    {
        elem = f(mesh[index]);
        ++index;
    }
    this->value = result;
}

double SliceFunction::operator[](size_t index) const
{
    return value[index];
}

SliceFunction operator+(const SliceFunction &f1, const SliceFunction &f2)
{
    auto value = f1.value;
    size_t index = 0;
    for (auto &elem : value)
    {
        elem += f2.value[index];
        ++index;
    }
    return SliceFunction(std::move(value));
}

SliceFunction operator-(const SliceFunction &f1, const SliceFunction &f2)
{
    auto value = f1.value;
    size_t index = 0;
    for (auto &elem : value)
    {
        elem -= f2.value[index];
        ++index;
    }
    return SliceFunction(std::move(value));
}

SliceFunction operator*(const SliceFunction &f1, const SliceFunction &f2)
{
    auto value = f1.value;
    size_t index = 0;
    for (auto &elem : value)
    {
        elem *= f2.value[index];
        ++index;
    }
    return SliceFunction(std::move(value));
}

SliceFunction operator/(const SliceFunction &f1, const SliceFunction &f2)
{
    auto value = f1.value;
    size_t index = 0;
    for (auto &elem : value)
    {
        elem /= f2.value[index];
        ++index;
    }
    return SliceFunction(std::move(value));
}

SliceFunction operator*(double scalar, const SliceFunction &f)
{
    auto value = f.value;
    size_t index = 0;
    for (auto &elem : value)
    {
        elem *= scalar;
        ++index;
    }
    return SliceFunction(std::move(value));
}

void SliceFunction::ApplyForeach(double (*function)(double))
{
    for (auto &elem : value)
    {
        elem = function(elem);
    }
}

Function::Function(const Mesh &mesh) : mesh(mesh)
{
}

void Function::SetInitialValue(double(* f_initial)(double), double initial_time)
{
    if (data.empty())
    {
        data.emplace_back(initial_time, SliceFunction(this->mesh, f_initial));
    }
}

void Function::SetInitialValue(SliceFunction &&initial_value, double initial_time)
{
    if (data.empty())
    {
        data.emplace_back(initial_time, initial_value);
    }
}

void Function::SetInitialValue(const SliceFunction &initial_value, double initial_time)
{
    if (data.empty())
    {
        data.emplace_back(initial_time, initial_value);
    }
}

double Function::GetLastTimeStep()
{
    return data.back().first;
}

SliceFunction Function::GetLastLayer()
{
    return data.back().second;
}

void Function::AddLayer(double time_step, SliceFunction &&value)
{
    double layer_time = this->GetLastTimeStep() + time_step;
    data.emplace_back(layer_time, value);
}