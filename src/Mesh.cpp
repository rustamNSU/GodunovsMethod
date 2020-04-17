#include "Mesh.hpp"

Mesh::Mesh() = default;
Mesh::Mesh(size_t points_number, double left, double right)
{
    this->SetMesh(points_number, left, right);
}

void Mesh::SetMesh(size_t points_number, double left, double right)
{
    this->mesh.resize(points_number);
    
    double step = (right - left) / (points_number - 1);
    for (auto &elem : this->mesh)
    {
        elem = left;
        left += step;
    }
    this->mesh.back() = right;
}

size_t Mesh::GetSize()
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
SliceFunction::SliceFunction(const Mesh &mesh)
{
    this->mesh = std::make_shared<Mesh>(mesh);
    this->value = std::vector<double>(this->mesh->GetSize(), 0.0);
}

SliceFunction::SliceFunction(std::shared_ptr<Mesh> mesh)
{
    this->mesh = mesh;
    this->value = std::vector<double>(mesh->GetSize(), 0.0);
}

SliceFunction::SliceFunction(const std::shared_ptr<Mesh> &mesh)
{
    this->mesh = mesh;
    this->value = std::vector<double>(mesh->GetSize(), 0.0);
}

void SliceFunction::SetValue(std::function<double(double)> f)
{
    int index = 0;
    for (auto &elem : this->value)
    {
        elem = f((*mesh)[index]);
        ++index;
    }
}

void Function::SetInitialValue(SliceFunction initial_value, double initial_time)
{
    if (data.size() != 0){
        return;
    }
    data.emplace_back(initial_time, initial_value);
}

void Function::SetInitialValue(const SliceFunction &initial_value, double initial_time)
{
    if (data.size() != 0){
        return;
    }
    data.push_back(std::pair(initial_time, initial_value));
}