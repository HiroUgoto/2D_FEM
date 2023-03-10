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
  input_wave::ricker(const Eigen::VectorXd tim, const double fp, const double tp, const double amp) {
    size_t num = tim.size();
    Eigen::VectorXd wave = Eigen::VectorXd::Zero(num);

    for (size_t i = 0 ; i < num ; i++ ){
      double t1 = std::pow((tim(i)-tp)*M_PI*fp,2.0);
      wave(i) = (2.0*t1-1.0) * std::exp(-t1) * amp;
    }

    return wave;
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

// ------------------------------------------------------------------- //
Eigen::VectorXd
  input_wave::tapered_sin(const Eigen::VectorXd tim, const double fp, const double taper, const double duration, const double amp) {
    size_t num = tim.size();
    Eigen::VectorXd wave = input_wave::simple_sin(tim,fp,amp);
    Eigen::VectorXd coeff = Eigen::VectorXd::Ones(num);

    for (size_t i = 0 ; i < num ; i++ ){
      if (tim(i) < taper) {
        coeff(i) = tim(i) / taper;
      } else if ((duration-taper < tim(i)) && (tim(i) <= duration)) {
        coeff(i) = (duration-tim(i)) / taper;
      } else if (duration < tim(i)) {
        coeff(i) = 0.0;
      }
    }

    return wave.array() * coeff.array();
  }

  // ------------------------------------------------------------------- //
Eigen::VectorXd
  input_wave::wavedata(const Eigen::VectorXd tim, const std::string inputfile) {
    std::ifstream f(inputfile);
    if (f.fail()) {
	    std::cout << "Cannot open file\n";
	    exit(0);
	  }

    std::string line;
    double number;
    size_t num = tim.size();
    Eigen::VectorXd wave = Eigen::VectorXd::Zero(num);
    for (size_t i = 0 ; i < num ; i++ ){
      std::getline(f,line);
      std::istringstream iss(line);
      iss >> number;
      wave(i) = number;
    }
    return wave;
  }