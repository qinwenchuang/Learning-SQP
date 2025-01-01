#include "../core/include/problem.hpp"

namespace example {

using namespace sqp;

class SimpleExample : public Problem {
 public:
  SimpleExample() = default;
  ~SimpleExample() = default;

  void Evaluate(const Vector &x, double *const value) {
    *value = log(1.0 + 100.0 * (x(1) - x(0) * x(0)) * (x(1) - x(0) * x(0)) +
                 (1.0 - x(0)) * (1.0 - x(0)));
  }

  void Jacobian(const Vector &x, Vector *const grad, double *const value) {
    const size_t N = x.rows();
    *grad = Vector::Zero(N);

    double temp = 1.0 + 100.0 * (x(1) - x(0) * x(0)) * (x(1) - x(0) * x(0)) +
                  (1.0 - x(0)) * (1.0 - x(0));
    (*grad)(0, 0) =
        (-400.0 * (x(1) - x(0) * x(0)) * x(0) - 2.0 * (1 - x(0))) / temp;
    (*grad)(1, 0) = (200.0 * (x(1) - x(0) * x(0))) / temp;

    *value = log(temp);
  }

  void Constraint(const Vector &x, Vector *const value, Vector *const l,
                  Vector *const u) {
    (*l) = Vector::Zero(1);
    (*u) = Vector::Zero(1);
    (*value)(0) = x(0) * x(0) + x(1) * x(1);

    (*l) << 0.0;
    (*u) << 1.0;
  }

  void ConstraintJacobian(const Vector &x, Matrix *const grad,
                          Vector *const value, Vector *const l,
                          Vector *const u) {
    (*l) = Vector::Zero(1);
    (*u) = Vector::Zero(1);

    (*grad)(0, 0) = 2.0 * x(0);
    (*grad)(1, 0) = 2.0 * x(1);
    (*value)(0) = x(0) * x(0) + x(1) * x(1);

    (*l) << 0.0;
    (*u) << 1.0;
  }
};

}  // namespace example