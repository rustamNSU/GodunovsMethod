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
    size_t GetSize() const;
    double operator[](size_t index) const;
    double& operator[](size_t index);
};


class SliceFunction
{
private:
    std::vector<double> value;

public:
    SliceFunction();
    SliceFunction(const Mesh &mesh, double(*f)(double));
    void SetValue(const Mesh &mesh, double(*f)(double));
};

class Function
{
private:
    /* Time + SliceFunction binding with time */
    std::vector<std::pair<double, SliceFunction>> data;
    const Mesh &mesh;

public:
    Function(const Mesh &mesh);

    void SetInitialValue(double(*f_initial)(double), double initial_time = 0.0);
    void SetInitialValue(SliceFunction &&initial_value, double initial_time = 0.0);
    void SetInitialValue(const SliceFunction &initial_value, double initial_time = 0.0);

};