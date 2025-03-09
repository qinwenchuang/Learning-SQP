#include <memory>

#include "core/include/integratior.hpp"
#include "core/include/sqp.hpp"
#include "example/simple_problem.hpp"

using namespace sqp;
using namespace example;

double test_sqp() {
  std::clock_t start_time, end_time;

  start_time = clock();

  size_t N = 2;
  size_t M = 1;
  size_t max_qp_iter = 100;
  size_t max_ls_iter = 100;
  Vector init_x = Vector::Zero(N);

  std::shared_ptr<Problem> simple_example = std::make_shared<SimpleExample>();

  std::shared_ptr<SQP> ptr_sqp =
      std::make_shared<SQP>(init_x, N, M, max_qp_iter, max_ls_iter);
  ptr_sqp->SetProblem(simple_example);
  ptr_sqp->Process();
  end_time = clock();

  std::cout << "X: " << ptr_sqp->GetResult()(0) << "  "
            << ptr_sqp->GetResult()(1) << std::endl;

  return (double)(end_time - start_time) / CLOCKS_PER_SEC;
}

void test_model_jac() {
  const size_t state_num = 3;
  const size_t control_num = 2;
  const double dt = 0.1;
  std::shared_ptr<Integrator> model_ptr =
      std::make_shared<Integrator>(state_num, control_num);
  Vector x = Vector::Zero(state_num);
  Vector u = Vector::Zero(control_num);
  Matrix out = Matrix::Zero(state_num, state_num + control_num);
  x << 0.0, 0.0, M_PI / 4.0;
  u << 10.0, 0.0;
  model_ptr->Jacobian(x, u, dt, &out);

  for (size_t i = 0; i < state_num; ++i) {
    std::cout << "row" + std::to_string(i) << ": " << out(i, 0) << " , "
              << out(i, 0) << " , " << out(i, 1) << " , " << out(i, 2) << " , "
              << out(i, 3) << " , " << out(i, 4) << std::endl;
  }
}

int main() {
  // test simple
  std::cout << " SQP: " << test_sqp() << std::endl;
  // test model jacobian
  test_model_jac();
  return 0;
}