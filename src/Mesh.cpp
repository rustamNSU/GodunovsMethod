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
    double element_center = left + step / 2.0;
    for (auto &elem : this->mesh)
    {
        elem = element_center;
        element_center += step;
    }
}

size_t Mesh::GetSize () const
{
    return mesh.size();
}

double Mesh::operator[](size_t index) const
{
    return mesh[index];
}

double& Mesh::operator[](size_t index)
{
    return mesh[index];
}

SliceFunction::SliceFunction() = default;

SliceFunction::SliceFunction(const Mesh &mesh, double(*f)(double))
{
    this->SetValue(mesh, f);
}

void SliceFunction::SetValue(const Mesh &mesh, double(*f)(double))
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

Function::Function(const Mesh &mesh): mesh(mesh){}

void Function::SetInitialValue(double(*f_initial)(double), double initial_time)
{
    if (data.empty()){
        data.emplace_back(initial_time, SliceFunction(this->mesh, f_initial));
    }
}
void Function::SetInitialValue(SliceFunction &&initial_value, double initial_time)
{
    if (data.empty()){
        data.emplace_back(initial_time, initial_value);
    }
}

void Function::SetInitialValue(const SliceFunction &initial_value, double initial_time)
{
    if (data.empty()){
        data.emplace_back(initial_time, initial_value);
    }
}