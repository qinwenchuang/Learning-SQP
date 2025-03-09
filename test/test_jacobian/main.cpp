#include <benchmark/benchmark.h>

#include <chrono>
#include <iostream>

#include "example_c.hpp"
#include "example_cpp.hpp"
#include "example_sym.hpp"

static void bench_casadi_c(benchmark::State& state) {
  double x[] = {0.0, 0.0, M_PI / 4.0};
  double u[] = {10.0, 0.0};
  double x1[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double x2[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double* args[] = {x, u};
  double* res[] = {x1, x2};
  casadi_int iw[400];
  casadi_real w[400];
  for (auto _ : state) {
    casadi_f0(args, res, iw, w, 0);
  }
}
BENCHMARK(bench_casadi_c);

static void bench_model_sym(benchmark::State& state) {
  Vector3 arr_x = {0.0, 0.0, M_PI / 4.0};
  Vector2 arr_u = {10.0, 0.0};
  Matrix33 jac_x;
  Matrix32 jac_u;
  TestSym test = TestSym();

  for (auto _ : state) {
    test.TestJacobianSym(arr_x, arr_u, jac_x, jac_u);
  }
}
BENCHMARK(bench_model_sym);

static void bench_model_cpp(benchmark::State& state) {
  Matrix cpp_out = Matrix::Zero(3, 5);
  Vector vec_x = Vector::Zero(3);
  Vector vec_u = Vector::Zero(2);
  vec_x << 0.0, 0.0, M_PI / 4.0;
  vec_u << 10.0, 0.0;
  double dt = 0.1;

  Integrator test = Integrator(3, 2);

  for (auto _ : state) {
    test.Jacobian(vec_x, vec_u, dt, &cpp_out);
  }
}
BENCHMARK(bench_model_cpp);

BENCHMARK_MAIN();

// int main() {
//   Vector3 arr_x = {0.0, 0.0, M_PI / 4.0};
//   Vector2 arr_u = {10.0, 0.0};
//   Matrix33 jac_x;
//   Matrix32 jac_u;
//   TestJacobianSym(arr_x, arr_u, jac_x, jac_u);
//   std::cout << jac_u[0][0] << " , " << jac_u[0][1] << " , " << jac_u[1][0]
//             << " , " << jac_u[1][1] << " , " << jac_u[2][0] << " , "
//             << jac_u[2][1] << std::endl;
// }