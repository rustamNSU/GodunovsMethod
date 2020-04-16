#pragma once

#include <vector>
#include <memory>

class Mesh
{
private:
    std::vector<double> mesh;

public:
    Mesh();
    Mesh(size_t points_number, double left, double right);
    void SetMesh(size_t points_number, double left, double right);
    size_t GetSize();
};


class DiscreteFunction
{
private:
    std::shared_ptr<Mesh> mesh;
    std::vector<double> value;

public:
    DiscreteFunction();
    DiscreteFunction(const Mesh &mesh);
    DiscreteFunction(const std::shared_ptr<Mesh> &mesh);

};