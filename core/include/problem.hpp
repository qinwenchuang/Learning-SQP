#pragma once

#include "common.hpp"

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

namespace sqp {

class Problem {
 public:
  Problem() = default;
  virtual ~Problem() = default;

  virtual void Evaluate(const Vector &x, double *value) {};
  virtual void Jacobian(const Vector &x, Vector *const grad,
                        double *const value) {};

  virtual void Constraint(const Vector &x, Vector *const value, Vector *const l,
                          Vector *const u) {};
  virtual void ConstraintJacobian(const Vector &x, Matrix *const grad,
                                  Vector *const value, Vector *const l,
                                  Vector *const u) {};
};

}  // namespace sqp