#include "all.h"
#include <Eigen/Core>
#include "material.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;

Material::Material () {}
Material::Material (size_t id, std::string style, std::vector<double> param) {
  this->id = id;
  this->style = style;
  this->param = param;

  this->set_param();
}

void Material::set_init(size_t id, std::string style, std::vector<double> param) {
  this->id = id;
  this->style = style;
  this->param = param;

  this->set_param();
}

void Material::print() {
  std::cout << this->id << ": ";
  std::cout << this->style << ", ";
  for (size_t i = 0 ; i < this->param.size() ; i++) {
    std::cout << this->param.at(i) << " ";
  }
  std::cout << "\n";
}

void Material::set_param() {
  if (this->style == "vs_vp_rho") {
    double vs = this->param.at(0);
    double vp = this->param.at(1);
    double rho = this->param.at(2);

    this->rmu = rho*vs*vs;
    this->rlambda = rho*vp*vp - 2.0*this->rmu;
    this->rho = rho;

  } else if (this->style == "nu_vp_rho") {
    double nu = this->param.at(0);
    double vp = this->param.at(1);
    double rho = this->param.at(2);

    this->rmu = rho/2.0*vp*vp*(1.0-2.0*nu)/(1.0-nu);
    this->rlambda = rho*nu*vp*vp/(1.0-nu);
    this->rho = rho;

  } else if (this->style == "nu_vs_rho") {
    double nu = this->param.at(0);
    double vs = this->param.at(1);
    double rho = this->param.at(2);

    this->rmu = rho*vs*vs;
    this->rlambda = 2.0*nu/(1.0-2.0*nu) * this->rmu;
    this->rho = rho;

  }
}

// ------------------------------------------------------------------- //
EM Material::mk_d(const size_t dof) {
    EM D;

    if (dof == 1) {
      D = EM::Zero(2,2);
      D(0,0) = this->rmu;
      D(1,1) = this->rmu;

    } else if (dof == 2) {
      D = EM::Zero(3,3);
      D(0,0) = this->rlambda + 2.0*this->rmu;
      D(0,1) = this->rlambda;
      D(1,0) = this->rlambda;
      D(1,1) = this->rlambda + 2.0*this->rmu;
      D(2,2) = this->rmu;

    } else if (dof == 2) {
      D = EM::Zero(5,5);
      D(0,0) = this->rlambda + 2.0*this->rmu;
      D(0,1) = this->rlambda;
      D(1,0) = this->rlambda;
      D(1,1) = this->rlambda + 2.0*this->rmu;
      D(2,2) = this->rmu;
      D(3,3) = this->rmu;
      D(4,4) = this->rmu;

    }

    return D;
  }

// ------------------------------------------------------------------- //
EM Material::mk_visco(const size_t dof) {
    EM D;
    double mu = 0.001; // [Pa s]

    if (dof == 1) {
      D = EM::Zero(2,2);
      D(0,0) = mu;
      D(1,1) = mu;

    } else if (dof == 2) {
      D = EM::Zero(3,3);
      D(2,2) = mu;

    } else if (dof == 2) {
      D = EM::Zero(5,5);
      D(2,2) = mu;
      D(3,3) = mu;
      D(4,4) = mu;

    }

    return D;
  }

// ------------------------------------------------------------------- //
EM Material::mk_imp(const size_t dof) {
    EM imp;
    double vs = sqrt(this->rmu/this->rho);
    double vp = sqrt((this->rlambda + 2.0*this->rmu)/this->rho);

    if (dof == 1) {
      imp = EM::Zero(1,1);
      imp(0,0) = this->rho * vs;

    } else if (dof == 2) {
      imp = EM::Zero(2,2);
      imp(0,0) = this->rho * vp;
      imp(1,1) = this->rho * vs;

    } else if (dof == 3) {
      imp = EM::Zero(3,3);
      imp(0,0) = this->rho * vp;
      imp(1,1) = this->rho * vs;
      imp(2,2) = this->rho * vs;
    }

    return imp;
  }
