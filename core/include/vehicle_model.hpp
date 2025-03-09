
#pragma once

#include "common.hpp"

namespace sqp {

class VehicleModel {
 public:
  VehicleModel() = default;
  virtual ~VehicleModel() = default;

  virtual void Evaluate(const Vector& x, const Vector& u, Vector& xdot) {
    double theta = x(2);
    double v = u(0);
    double delta = u(1);
    xdot(0) = v * cos(theta);
    xdot(1) = v * sin(theta);
    xdot(2) = v / L * tan(delta);
  }

  virtual void Jacobian(const Vector& x, const Vector& u, Matrix& state_jac,
                        Matrix& control_jac) {
    double theta = x(2);
    double v = u(0);
    double delta = u(1);
    state_jac(0, 2) = -v * sin(theta);
    state_jac(1, 2) = v * cos(theta);

    control_jac(0, 0) = cos(theta);
    control_jac(1, 0) = sin(theta);
    control_jac(2, 0) = 1.0 / L * tan(delta);
    control_jac(2, 1) = v / (L * cos(delta) * cos(delta));

  }

 private:
  double L = 2.4;
};

}  // namespace sqp