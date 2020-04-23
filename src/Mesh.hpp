#pragma once

#include <vector>
#include <memory>
#include <utility>

class Mesh {
private:
    std::vector<double> mesh;
    double step_x; // The delta x

public:
    Mesh();

    Mesh(size_t points_number, double left, double right);

    void SetMesh(size_t points_number, double left, double right);

    size_t GetSize() const;

    double GetStep() const;

    double operator[](size_t index) const;

    double &operator[](size_t index);
};


class SliceFunction {
private:
    std::vector<double> value;

public:
    SliceFunction();

    explicit SliceFunction(std::vector<double> &&value);

    SliceFunction(const Mesh &mesh, double(* f)(double));

    void SetValue(const Mesh &mesh, double(* f)(double));

    double operator[](size_t index) const;

    friend SliceFunction operator+(const SliceFunction &f1, const SliceFunction &f2);

    friend SliceFunction operator-(const SliceFunction &f1, const SliceFunction &f2);

    friend SliceFunction operator*(const SliceFunction &f1, const SliceFunction &f2);

    friend SliceFunction operator/(const SliceFunction &f1, const SliceFunction &f2);

    friend SliceFunction operator*(double scalar, const SliceFunction &f);

    void ApplyForeach(double (*function)(double));
};

class Function {
private:
    /* Time + SliceFunction binding with time */
    std::vector<std::pair<double, SliceFunction>> data;
    const Mesh &mesh;

public:
    Function(const Mesh &mesh);

    void SetInitialValue(double(* f_initial)(double), double initial_time = 0.0);

    void SetInitialValue(SliceFunction &&initial_value, double initial_time = 0.0);

    void SetInitialValue(const SliceFunction &initial_value, double initial_time = 0.0);

    double GetLastTimeStep();

    SliceFunction GetLastLayer();

    void AddLayer(double time_step, SliceFunction &&value);

};