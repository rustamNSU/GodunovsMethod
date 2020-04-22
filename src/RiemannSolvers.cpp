#include "RiemannSolvers.hpp"

EulerVariables::EulerVariables(
        double center_coordinate,
        double velocity,
        double density,
        double pressure,
        double gamma) {
    this->u = velocity;
    this->rho = density;
    this->pressure = pressure;
    this->x = center_coordinate;
    this->gamma = gamma;
    this->sound_speed = std::sqrt(gamma * pressure / rho);
    this->enthalpy = u * u / 2.0 + sound_speed * sound_speed / (gamma - 1.0);
    this->energy = rho * enthalpy - pressure;
}

VariableVector3 EulerVariables::ConservativeVariable() {
    return VariableVector3(
            x,
            rho,
            rho * u,
            energy
    );
}

VariableVector3 operator+(const VariableVector3 &v1, const VariableVector3 &v2) {
    VariableVector3 result(
            (v1.x + v2.x) / 2.0,
            v1.q1 + v2.q1,
            v1.q2 + v2.q2,
            v1.q3 + v2.q3
    );
    return result;
}

VariableVector3 operator-(const VariableVector3 &v1, const VariableVector3 &v2) {
    VariableVector3 result(
            (v1.x - v2.x) / 2.0,
            v1.q1 - v2.q1,
            v1.q2 - v2.q2,
            v1.q3 - v2.q3
    );
    return result;
}

VariableVector3 operator*(double scalar, const VariableVector3 &v) {
    VariableVector3 result(
            v.x,
            scalar * v.q1,
            scalar * v.q2,
            scalar * v.q3
    );
    return result;
}

VariableVector3 CalculateFlux(
        double velocity,
        double enthalpy,
        double gamma,
        const VariableVector3 &conservative_variable) {
    VariableVector3 flux;
    flux.x = conservative_variable.x;
    flux.q1 = conservative_variable.q2;
    flux.q2 = (gamma - 3.0) * velocity * velocity / 2.0 * conservative_variable.q1 +
              (3.0 - gamma) * velocity * conservative_variable.q2 +
              (gamma - 1) * conservative_variable.q3;
    flux.q3 = ((gamma - 1.0) * velocity * velocity * velocity / 2.0 - velocity * enthalpy) * conservative_variable.q1 +
              (enthalpy - (gamma - 1) * velocity * velocity) * conservative_variable.q2 +
              gamma * velocity * conservative_variable.q3;
    return flux;
}

VariableVector3 RoeSolver(EulerVariables f_left, EulerVariables f_right) {
    double sqrt_rho_left = std::sqrt(f_left.rho);
    double sqrt_rho_right = std::sqrt(f_right.rho);

    /* Roe averages */
    double average_velocity = (sqrt_rho_left * f_left.u + sqrt_rho_right * f_right.u) /
                              (sqrt_rho_left + sqrt_rho_right);
    double average_enthalpy = (sqrt_rho_left * f_left.enthalpy + sqrt_rho_right * f_right.enthalpy) /
                              (sqrt_rho_left + sqrt_rho_right);
    double average_sound_speed = std::sqrt(
            (f_left.gamma - 1) *
            (average_enthalpy - average_velocity * average_velocity / 2.0));
    double average_rho = sqrt_rho_left * sqrt_rho_right;

    double delta_pressure = f_right.pressure - f_left.pressure;
    double delta_velocity = f_right.u - f_left.u;
    double aux_coefficient = 1.0 / average_sound_speed / average_sound_speed;

    double alpha1 = (delta_pressure - average_sound_speed * average_rho * delta_velocity) * aux_coefficient / 2.0;
    double alpha2 = (f_right.rho - f_left.rho) - aux_coefficient * delta_pressure;
    double alpha3 = (delta_pressure + average_sound_speed * average_rho * delta_velocity) * aux_coefficient / 2.0;

    double eigen_value1 = average_velocity - average_sound_speed;
    double eigen_value2 = average_velocity;
    double eigen_value3 = average_velocity + average_sound_speed;

    double x_boundary = (f_left.x + f_right.x) / 2.0;
    VariableVector3 r1(
            x_boundary,
            1.0,
            average_velocity - average_sound_speed,
            average_enthalpy - average_velocity * average_sound_speed
    );
    VariableVector3 r2(
            x_boundary,
            1.0,
            average_velocity,
            average_velocity * average_velocity / 2.0
    );
    VariableVector3 r3(
            x_boundary,
            1.0,
            average_velocity + average_sound_speed,
            average_enthalpy + average_velocity * average_sound_speed
    );

    VariableVector3 flux;
    flux = 0.5 * CalculateFlux(average_velocity, average_enthalpy, f_left.gamma, f_left.ConservativeVariable() +
                                                                                 f_right.ConservativeVariable()) -
           0.5 * (std::abs(eigen_value1) * alpha1 * r1 +
                  std::abs(eigen_value2) * alpha2 * r2 +
                  std::abs(eigen_value3) * alpha3 * r3);

    return flux;
}