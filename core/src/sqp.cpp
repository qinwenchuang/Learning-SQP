#include "sqp.hpp"

#include <qpOASES/QProblem.hpp>
#include <qpOASES/QProblemB.hpp>

namespace sqp {

SQP::SQP(Vector x, size_t n, const size_t m, const size_t qp_iter,
         const size_t ls_iter)
    : N_(n), M_(m), max_qp_iter_(qp_iter), max_ls_iter_(ls_iter), x_(x) {
  direction_ = Vector::Zero(N_);
  qpoases_iter_ = 0;
  ls_iter_ = 0;
  sqp_iter_ = 0;

  step_ = Vector::Zero(N_);
  grad_ = Vector::Zero(N_);
  obj_ = 0.0;
  con_ = Vector::Zero(M_);
  l_ = Vector::Zero(M_);
  u_ = Vector::Zero(M_);
  grad_obj_ = Vector::Zero(N_);
  grad_con_ = Matrix::Zero(N_, N_);
  delta_grad_ = Vector::Zero(N_);
  Hessian_ = Matrix::Identity(N_, N_);
}

bool SQP::Process() {
  bool flag = true;

  size_t i = 0;
  for (i = 0; i < max_qp_iter_; ++i) {
    // direction
    CalculateQP(i);
    // linesearch
    double alpha = CalculateLS(i);
    // take step
    x_ = x_ + alpha * direction_;
    step_ = alpha * direction_;

    if (abs(step_.maxCoeff()) < 1e-7 && i > 1) {
      break;
    }

    if (abs(delta_grad_.maxCoeff()) < 1e-7 && i > 1) {
      break;
    }
  }
  sqp_iter_ = i;
  return flag;
}

bool SQP::CalculateQP(const size_t& iter) {
  problem_ptr_->Jacobian(x_, &grad_obj_, &obj_);
  problem_ptr_->ConstraintJacobian(x_, &grad_con_, &con_, &l_, &u_);

  if (iter < 1) {
    Hessian_ = Matrix::Identity(N_, N_);
  } else {
    Vector new_grad = grad_obj_;
    delta_grad_ = new_grad - grad_;
    Common::LBFGS(Hessian_, step_, delta_grad_);
  }
  grad_ = grad_obj_;

  double H[N_ * N_];
  double g[N_];
  double A[M_ * N_];
  double lb[N_];
  double ub[N_];
  double lbA[M_];
  double ubA[M_];

  for (size_t i = 0; i < N_; ++i) {
    for (size_t j = 0; j < N_; ++j) {
      H[i * N_ + j] = Hessian_(i, j);
    }
    g[i] = grad_(i);
    lb[i] = -std::numeric_limits<double>::max();  // Note: Adapt simple exampe
    ub[i] = std::numeric_limits<double>::max();   // Note: Adapt simple exampe
  }

  for (size_t i = 0; i < M_; ++i) {
    for (size_t j = 0; j < N_; ++j) {
      A[i * N_ + j] = grad_con_(i, j);
    }
    lbA[i] = l_(i) - con_(i);
    ubA[i] = u_(i) - con_(i);
  }

  qpOASES::QProblem qp_solver(N_, M_);
  qpOASES::Options options;
  options.enableCholeskyRefactorisation = 1;
  options.setToMPC();
  options.printLevel = qpOASES::PrintLevel::PL_NONE;
  qp_solver.setOptions(options);

  int qp_iter = 100;
  qp_solver.init(H, g, A, lb, ub, lbA, ubA, qp_iter);
  qpoases_iter_ += qp_iter + 1;

  double xOpt[N_];
  qp_solver.getPrimalSolution(xOpt);
  for (size_t i = 0; i < N_; ++i) {
    direction_(i) = xOpt[i];
  }
  return true;
}

double SQP::CalculateLS(const size_t& iter) {
  bool flag = true;
  double rho = 0.5;
  double alpha = 1.00;
  double eta = 0.25;
  double tau = 0.5;

  double norm = ConsNorm(con_, l_, u_);

  double mu = (grad_obj_.dot(direction_) +
               0.5 * direction_.dot(Hessian_ * direction_)) /
              ((1.0 - rho) * norm);

  double phi = obj_ + mu * norm;
  double dp_phi = grad_obj_.dot(direction_) - mu * norm;

  size_t i = 0;
  for (i = 0; i < max_ls_iter_; ++i) {
    double temp_obj = 0.0;
    Vector temp_x = x_ + alpha * direction_;

    problem_ptr_->Evaluate(temp_x, &temp_obj);
    problem_ptr_->Constraint(temp_x, &con_, &l_, &u_);

    double temp_phi = temp_obj + mu * ConsNorm(con_, l_, u_);

    if (temp_phi <= phi + alpha * eta * dp_phi) {
      break;
    } else {
      alpha *= tau;
    }
    ls_iter_ += i + 1;
  }
  return alpha;
}

double SQP::ConsNorm(const Vector& cons, const Vector& l, const Vector& u) {
  double norm = (l - cons).cwiseAbs().sum() + (cons - u).cwiseAbs().sum();
  return norm;
}

}  // namespace sqp