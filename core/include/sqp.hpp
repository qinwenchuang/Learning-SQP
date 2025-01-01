#pragma once

#include "common.hpp"
#include "problem.hpp"
#include <limits>

namespace sqp {

class SQP {
 public:
  SQP(Vector x, size_t n, const size_t m, const size_t qp_iter,
      const size_t ls_iter);
  ~SQP() = default;

  bool Process();
  bool CalculateQP(const size_t& iter);
  double CalculateLS(const size_t& iter);

  const Vector& GetResult() const { return x_; };
  size_t GetQpIter() const { return qpoases_iter_; };
  size_t GetLsIter() const { return ls_iter_; };
  size_t GetSqpIter() const { return sqp_iter_; };
  double GetStepNorm() const { return step_.norm(); };

  void SetProblem(std::shared_ptr<Problem> ptr) { problem_ptr_ = ptr; };

 private:
  double ConsNorm(const Vector& cons, const Vector& l, const Vector& u);

 private:
  size_t N_ = 0;
  size_t M_ = 0;
  Vector x_;

  size_t max_qp_iter_ = 0;
  size_t max_ls_iter_ = 0;

  Vector step_;
  Vector grad_;

  double obj_ = 0.0;
  Vector con_;
  Vector l_;
  Vector u_;
  Vector grad_obj_;
  Matrix grad_con_;

  Vector delta_grad_;
  Matrix Hessian_;
  Vector direction_;

  size_t qpoases_iter_ = 0;
  size_t ls_iter_ = 0;
  size_t sqp_iter_ = 0;

  std::shared_ptr<Problem> problem_ptr_ = std::make_shared<Problem>();
};

}  // namespace sqp
