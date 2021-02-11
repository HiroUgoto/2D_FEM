#include "all.h"
#include <Eigen/Core>
#include "material.h"

Material::Material (size_t id, std::string style, std::vector<double> param) {
    Material::id = id;
    Material::style = style;
    Material::param = param;

    Material::set_param();
}

void Material::print() {
  std::cout << Material::id << ": ";
  std::cout << Material::style << ", ";
  for (size_t i = 0 ; i < Material::param.size() ; i++) {
    std::cout << Material::param.at(i) << " ";
  }
  std::cout << "\n";
}

void Material::set_param() {
  if (Material::style == "vs_vp_rho") {
    double vs = Material::param.at(0);
    double vp = Material::param.at(1);
    double rho = Material::param.at(2);

    Material::rmu = rho*vs*vs;
    Material::rlambda = rho*vp*vp - 2.0*Material::rmu;
    Material::rho = rho;

  } else if (Material::style == "nu_vp_rho") {
    double nu = Material::param.at(0);
    double vp = Material::param.at(1);
    double rho = Material::param.at(2);

    Material::rmu = rho/2.0*vp*vp*(1.0-2.0*nu)/(1.0-nu);
    Material::rlambda = rho*nu*vp*vp/(1.0-nu);
    Material::rho = rho;

  } else if (Material::style == "nu_vs_rho") {
    double nu = Material::param.at(0);
    double vs = Material::param.at(1);
    double rho = Material::param.at(2);

    Material::rmu = rho*vs*vs;
    Material::rlambda = 2.0*nu/(1.0-2.0*nu) * Material::rmu;
    Material::rho = rho;

  }
}

// ------------------------------------------------------------------- //
Eigen::MatrixXd
  Material::mk_d(const size_t dof) {
  Eigen::MatrixXd D;

  if (dof == 1) {
    D = Eigen::MatrixXd::Zero(2,2);
    D(0,0) = Material::rmu;
    D(1,1) = Material::rmu;

  } else if (dof == 2) {
    D = Eigen::MatrixXd::Zero(3,3);
    D(0,0) = Material::rlambda + 2.0*Material::rmu;
    D(0,1) = Material::rlambda;
    D(1,0) = Material::rlambda;
    D(1,1) = Material::rlambda + 2.0*Material::rmu;
    D(2,2) = Material::rmu;

  } else if (dof == 2) {
    D = Eigen::MatrixXd::Zero(5,5);
    D(0,0) = Material::rlambda + 2.0*Material::rmu;
    D(0,1) = Material::rlambda;
    D(1,0) = Material::rlambda;
    D(1,1) = Material::rlambda + 2.0*Material::rmu;
    D(2,2) = Material::rmu;
    D(3,3) = Material::rmu;
    D(4,4) = Material::rmu;

  }

  return D;
}

// ------------------------------------------------------------------- //
Eigen::MatrixXd
  Material::mk_visco(const size_t dof) {
  Eigen::MatrixXd D;
  double mu = 0.001; // [Pa s]

  if (dof == 1) {
    D = Eigen::MatrixXd::Zero(2,2);
    D(0,0) = mu;
    D(1,1) = mu;

  } else if (dof == 2) {
    D = Eigen::MatrixXd::Zero(3,3);
    D(2,2) = mu;

  } else if (dof == 2) {
    D = Eigen::MatrixXd::Zero(5,5);
    D(2,2) = mu;
    D(3,3) = mu;
    D(4,4) = mu;

  }

  return D;
}
