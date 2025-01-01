#include <memory>

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

int main() {
  std::cout << " SQP: " << test_sqp() << std::endl;
  return 0;
}