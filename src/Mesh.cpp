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

DiscreteFunction::DiscreteFunction() = default;
DiscreteFunction::DiscreteFunction(const Mesh &mesh)
{
    this->mesh = std::make_shared<Mesh>(mesh);
    this->value = std::vector<double>(this->mesh->GetSize(), 0.0);
}

DiscreteFunction::DiscreteFunction(const std::shared_ptr<Mesh> &mesh)
{
    this->mesh = mesh;
    this->value = std::vector<double>(mesh->GetSize(), 0.0);
}