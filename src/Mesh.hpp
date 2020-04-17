#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <utility>

class Mesh
{
private:
    std::vector<double> mesh;

public:
    Mesh();
    Mesh(size_t points_number, double left, double right);
    void SetMesh(size_t points_number, double left, double right);
    size_t GetSize();
    double operator[](size_t index) const;
    double& operator[](size_t index);
};


class SliceFunction
{
private:
    std::shared_ptr<Mesh> mesh;
    std::vector<double> value;

public:
    SliceFunction();
    SliceFunction(const Mesh &mesh);
    SliceFunction(std::shared_ptr<Mesh> mesh);
    SliceFunction(const std::shared_ptr<Mesh> &mesh);

    void SetValue(std::function<double(double)> f);
};

class Function
{
private:
    /* Time + SliceFunction bbinding with time */
    std::vector<std::pair<double, SliceFunction>> data;
    std::shared_ptr<Mesh> mesh;

public:
    Function();

    void SetInitialValue(SliceFunction initial_value, double initial_time = 0.0);
    void SetInitialValue(const SliceFunction &initial_value, double initial_time = 0.0);

};