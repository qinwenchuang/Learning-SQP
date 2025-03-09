#pragma once

#include "common.hpp"
#include "problem.hpp"
#include "vehicle_model.hpp"

namespace sqp {
class Integrator {
 public:
  Integrator(size_t n, size_t m) : N_(n), M_(m) {};
  ~Integrator() = default;

  virtual void Evaluate(const Vector &x, const Vector &u, const double T,
                        Vector *value);
  virtual void Jacobian(const Vector &x, const Vector &u, const double T,
                        Matrix *value);

 private:
  size_t N_ = 0;
  size_t M_ = 0;

  Vector K1_ = Vector::Zero(N_);
  Vector K2_ = Vector::Zero(N_);
  Vector K3_ = Vector::Zero(N_);
  Vector K4_ = Vector::Zero(N_);

  Matrix A0_ = Matrix::Zero(N_, N_);
  Matrix B0_ = Matrix::Zero(N_, M_);
  Matrix A1_ = Matrix::Zero(N_, N_);
  Matrix B1_ = Matrix::Zero(N_, M_);
  Matrix A2_ = Matrix::Zero(N_, N_);
  Matrix B2_ = Matrix::Zero(N_, M_);
  Matrix A3_ = Matrix::Zero(N_, N_);
  Matrix B3_ = Matrix::Zero(N_, M_);
  std::shared_ptr<VehicleModel> model_ptr = std::make_shared<VehicleModel>();
};

}  // namespace sqp