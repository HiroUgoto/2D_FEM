#include "all.h"
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include "elasto_plastic.h"
#include "ep_model.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;


// ------------------------------------------------------------------- //
EP* set_ep_style(double dof, std::string style, std::vector<double> param_ep) {
  EP* ep_p = nullptr;

  if (style == "ep_Li") {
    double nu = param_ep.at(0);
    double G0 = param_ep.at(1);
    double M  = param_ep.at(2);
    double e0 = param_ep.at(3);
    double eg = param_ep.at(4);
    double d1 = param_ep.at(5);

    ep_p = new Li2002(G0,nu,M,eg,d1);
    ep_p->e0 = e0;

    ep_p->dof = dof;
    ep_p->style = style;
    ep_p->clear_strain();
  }

  return ep_p;
}

// ------------------------------------------------------------------- //
EP::EP () {}
EP::~EP () {}

// ------------------------------------------------------------------- //
EM EP::FEMstress_to_matrix(EV FEMstress) {
  EV stress_vec = EV::Zero(6);
  if (this->dof == 1) {
    stress_vec[0] = 0.0         ; stress_vec[1] = 0.0;
    stress_vec[2] = 0.0         ; stress_vec[3] = FEMstress[0];
    stress_vec[4] = FEMstress[1]; stress_vec[5] = 0.0;
  } else if (this->dof == 2) {
    stress_vec[0] = -FEMstress[0]; stress_vec[1] = -FEMstress[0];
    stress_vec[2] = -FEMstress[1]; stress_vec[3] = 0.0;
    stress_vec[4] =  0.0         ; stress_vec[5] = FEMstress[2];
  } else if (this->dof == 3) {
    stress_vec[0] = -FEMstress[0]; stress_vec[1] = -FEMstress[0];
    stress_vec[2] = -FEMstress[1]; stress_vec[3] =  FEMstress[3];
    stress_vec[4] =  FEMstress[4]; stress_vec[5] =  FEMstress[2];
  }

  return this->vector_to_matrix(stress_vec);
}

EM EP::FEMstrain_to_matrix(EV FEMstrain) {
  EV strain_vec = EV::Zero(6);
  if (this->dof == 1) {
    strain_vec[0] = 0.0         ; strain_vec[1] = 0.0;
    strain_vec[2] = 0.0         ; strain_vec[3] = 0.5*FEMstrain[0];
    strain_vec[4] = 0.5*FEMstrain[1]; strain_vec[5] = 0.0;
  } else if (this->dof == 2) {
    strain_vec[0] = -FEMstrain[0]; strain_vec[1] = 0.0;
    strain_vec[2] = -FEMstrain[1]; strain_vec[3] = 0.0;
    strain_vec[4] =  0.0         ; strain_vec[5] = 0.5*FEMstrain[2];
  } else if (this->dof == 3) {
    strain_vec[0] = -FEMstrain[0]; strain_vec[1] = 0.0;
    strain_vec[2] = -FEMstrain[1]; strain_vec[3] = 0.5*FEMstrain[3];
    strain_vec[4] = 0.5*FEMstrain[4]; strain_vec[5] = 0.5*FEMstrain[2];
  }

  return this->vector_to_matrix(strain_vec);
}

// ------------------------------------------------------------------- //
EV EP::matrix_to_FEMstress(EM stress) {
  EV FEMstress;
  EV stress_vec = this->matrix_to_vector(stress);
  if (this->dof == 1) {
    FEMstress = EV::Zero(2);
    FEMstress(0) = stress_vec(3); FEMstress(1) = stress_vec(4);
  } else if (this->dof == 2) {
    FEMstress = EV::Zero(3);
    FEMstress(0) = -stress_vec(0); FEMstress(1) = -stress_vec(2);
    FEMstress(2) =  stress_vec(5);
  } else if (this->dof == 3) {
    FEMstress = EV::Zero(5);
    FEMstress(0) = -stress_vec(0); FEMstress(1) = -stress_vec(2);
    FEMstress(2) = stress_vec(5);
    FEMstress(3) = stress_vec(3); FEMstress(4) = stress_vec(4);
  }
  return FEMstress;
}

// ------------------------------------------------------------------- //
EM EP::modulus_to_Dmatrix(Eigen::Tensor<double,4> E) {
  EM D;
  if (this->dof == 1) {
    D = EM::Zero(2,2);
    D(0,0) = E(0,1,0,1); D(0,1) = E(0,1,1,2);
    D(1,0) = E(1,2,0,1); D(1,1) = E(1,2,1,2);
  } else if (this->dof == 2) {
    D = EM::Zero(3,3);
    D(0,0) = E(0,0,0,0); D(0,1) = E(0,0,2,2); D(0,2) = 0.5*(E(0,0,2,0)+E(0,0,0,2));
    D(1,0) = E(2,2,0,0); D(1,1) = E(2,2,2,2); D(1,2) = 0.5*(E(2,2,2,0)+E(2,2,0,2));

    D(2,0) = 0.5*(E(2,0,0,0)+E(0,2,0,0));
    D(2,1) = 0.5*(E(2,0,2,2)+E(0,2,2,2));
    D(2,2) = 0.25*(E(2,0,2,0)+E(0,2,2,0)+E(2,0,0,2)+E(0,2,0,2));

    // D(0,0) = E(0,0,0,0); D(0,1) = E(0,0,2,2); D(0,2) = E(0,0,2,0);
    // D(1,0) = E(2,2,0,0); D(1,1) = E(2,2,2,2); D(1,2) = E(2,2,2,0);
    // D(2,0) = E(2,0,0,0);
    // D(2,1) = E(2,0,2,2);
    // D(2,2) = E(2,0,2,0);
  }
  return D;
}

// ------------------------------------------------------------------- //
EM EP::vector_to_matrix(EV vec) {
  EM mat = EM::Zero(3,3);

  mat(0,0) = vec(0); mat(0,1) = vec(3); mat(0,2) = vec(5);
  mat(1,0) = vec(3); mat(1,1) = vec(1); mat(1,2) = vec(4);
  mat(2,0) = vec(5); mat(2,1) = vec(4); mat(2,2) = vec(2);
  return mat;
}

EV EP::matrix_to_vector(EM mat) {
  EV vec = EV::Zero(6);

  vec(0) = mat(0,0); vec(1) = mat(1,1); vec(2) = mat(2,2);
  vec(3) = 0.5*(mat(0,1)+mat(1,0));
  vec(4) = 0.5*(mat(1,2)+mat(2,1));
  vec(5) = 0.5*(mat(2,0)+mat(0,2));

  return vec;
}

// ------------------------------------------------------------------- //
