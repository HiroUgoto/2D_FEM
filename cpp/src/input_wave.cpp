#include "all.h"
#include <Eigen/Core>
#include "input_wave.h"

// ------------------------------------------------------------------- //
std::tuple<Eigen::VectorXd, double>
  input_wave::linspace(const double start, const double end, const size_t num) {
    double d = (end - start)/num;
    Eigen::VectorXd array = Eigen::VectorXd::Zero(num);

    for (size_t i = 0 ; i < num ; i++ ){
      array(i) = start + d*i;
    }

    return std::forward_as_tuple(array,d);
  }

// ------------------------------------------------------------------- //
Eigen::VectorXd
  input_wave::simple_sin(const Eigen::VectorXd tim, const double fp, const double amp) {
    size_t num = tim.size();
    Eigen::VectorXd wave = Eigen::VectorXd::Zero(num);

    for (size_t i = 0 ; i < num ; i++ ){
      wave(i) = amp * std::sin(2.0*fp * tim(i) * M_PI);
    }

    return wave;
  }
