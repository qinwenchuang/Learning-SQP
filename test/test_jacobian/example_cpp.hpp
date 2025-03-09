#include <Eigen/Dense>
#include <cmath>
#include <memory>

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

class VehicleModel {
 public:
  VehicleModel() = default;
  virtual ~VehicleModel() = default;

  virtual void Evaluate(const Vector &x, const Vector &u, Vector &xdot) {
    double theta = x(2);
    double v = u(0);
    double delta = u(1);
    xdot(0) = v * cos(theta);
    xdot(1) = v * sin(theta);
    xdot(2) = v / L * tan(delta);
  }

  virtual void Jacobian(const Vector &x, const Vector &u, Matrix &state_jac,
                        Matrix &control_jac) {
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

class Integrator {
 public:
  Integrator(size_t n, size_t m) : N_(n), M_(m) {};
  ~Integrator() = default;

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