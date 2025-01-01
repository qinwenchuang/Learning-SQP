#pragma once

#include <Eigen/Dense>
#include <ctime>
#include <iostream>
#include <limits>
#include <memory>
#include <limits>

#include "problem.hpp"

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

namespace sqp {

class Common {
 public:
  Common() = default;
  ~Common() = default;

  static void LBFGS(Matrix& B, const Vector& s, const Vector& y) {
    Vector Bs = B * s;
    Vector r = y;
    double sBs = s.dot(Bs);
    double sy = s.dot(y);
    double sr = sy;

    if (sy < 0.2 * sBs) {
      double theta = 0.8 * sBs / (sBs - sy);
      r = theta * y + (1.0 - theta) * Bs;
      sr = theta * sy + (1.0 - theta) * sBs;
    } else {
      r = y;
      sr = sy;
    }

    if (sr < std::numeric_limits<double>::epsilon()) {
      return;
    }

    B.noalias() += -(Bs * Bs.transpose()) / sBs + r * r.transpose() / sr;
  };
};

}  // namespace sqp
