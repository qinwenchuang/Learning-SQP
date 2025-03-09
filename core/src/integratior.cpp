#include "integratior.hpp"

namespace sqp {

void Integrator::Evaluate(const Vector &x, const Vector &u, const double T,
                          Vector *value) {
  model_ptr->Evaluate(x, u, K1_);
  model_ptr->Evaluate(x + K1_ * 0.5 * T, u, K2_);
  model_ptr->Evaluate(x + K2_ * 0.5 * T, u, K3_);
  model_ptr->Evaluate(x + K3_ * T, u, K4_);
  *value = x + T * (K1_ + 2 * K2_ + 2 * K3_ + K4_) / 6.0;
}

void Integrator::Jacobian(const Vector &x, const Vector &u, const double T,
                          Matrix *value) {
  model_ptr->Evaluate(x, u, K1_);
  model_ptr->Evaluate(x + K1_ * 0.5 * T, u, K2_);
  model_ptr->Evaluate(x + K2_ * 0.5 * T, u, K3_);

  model_ptr->Jacobian(x, u, A0_, B0_);
  model_ptr->Jacobian(x + 0.5 * K1_ * T, u, A1_, B1_);
  model_ptr->Jacobian(x + 0.5 * K2_ * T, u, A2_, B2_);
  model_ptr->Jacobian(x + K3_ * T, u, A3_, B3_);

  Matrix A0T = A0_ * T;
  Matrix A1T = A1_ * (Matrix::Identity(N_, N_) + 0.5 * A0T) * T;
  Matrix A2T = A2_ * (Matrix::Identity(N_, N_) + 0.5 * A1T) * T;
  Matrix A3T = A3_ * (Matrix::Identity(N_, N_) + A2T) * T;

  Matrix B0T = B0_ * T;
  Matrix B1T = B1_ * T + 0.5 * A1_ * B0T * T;
  Matrix B2T = B2_ * T + 0.5 * A2_ * B1T * T;
  Matrix B3T = B3_ * T + A3_ * B2T * T;

  value->block(0, 0, N_, N_) =
      Matrix::Identity(N_, N_) + (A0T + 2.0 * A1T + 2.0 * A2T + A3T) / 6.0;
  value->block(0, N_, N_, M_) = (B0T + 2.0 * B1T + 2.0 * B2T + B3T) / 6.0;
}

}  // namespace sqp