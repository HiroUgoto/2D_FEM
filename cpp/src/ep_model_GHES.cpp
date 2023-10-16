#include "all.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>
#include "elasto_plastic.h"
#include "ep_model_GHES.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;

GHE_S::GHE_S (double G0, double gr, double hmax, double nu) {
    // parameters
    this->G0 = G0;
    this->gr = gr;
    this->hmax = hmax;
    this->nu = nu;

    // stress and strain
    this->stress = EM::Zero(3,3);
    this->strain = EM::Zero(3,3);

}

// ------------------------------------------------------------------- //
void GHE_S::print() {
  // std::cout << this->G0 << "\n";
}

void GHE_S::clear_strain() {
  this->strain = EM::Zero(3,3);
}

// ------------------------------------------------------------------- //
std::tuple<double, double>
  GHE_S::elastic_modulus(const double e, const double p) {
    double G = this->G0;
    double K = G*2*(1+this->nu)/(3*(1-2*this->nu));
    return {G, K};
  }

std::tuple<double, double>
  GHE_S::elastic_modulus_lame() {
    auto [G0, K0] = this->elastic_modulus(0.0,0.0);
    return {G0, K0 - G0*2.0/3.0};
  }


// ------------------------------------------------------------------- //
void GHE_S::initial_state_isotropic(EV init_stress) {
//   // isotropic compression
//   double compression_stress = std::max(-init_stress(0),-init_stress(1));
//   this->isotropic_compression(this->e0,compression_stress);
//   this->e = this->e0;

//   // set initial parameters
//   double p = this->set_stress_variable_p(this->stress);
//   this->beta = p;
//   this->H2 = p;

}

void GHE_S::initial_state(EV init_stress) {
//   EM init_stress_mat = this->FEMstress_to_matrix(init_stress);

//   // isotropic compression
//   double compression_stress = std::max(-init_stress(0),-init_stress(1));
//   this->isotropic_compression(this->e0,compression_stress);
//   this->e = this->e0;

//   // set initial parameters
//   double p = this->set_stress_variable_p(this->stress);
//   this->beta = p;
//   this->H2 = p;

//   size_t nstep = 50;
//   EV dstrain_vec = EV::Zero(6);
//   EV dstress_vec = EV::Zero(6);

//   dstress_vec(0) = (init_stress_mat(0,0)-compression_stress)/nstep;
//   dstress_vec(1) = (init_stress_mat(1,1)-compression_stress)/nstep;
//   dstress_vec(2) = (init_stress_mat(2,2)-compression_stress)/nstep;

//   EM dstrain_input = this->vector_to_matrix(dstrain_vec);
//   EM dstress_input = this->vector_to_matrix(dstress_vec);

//   StateParameters sp0(this->strain,this->stress,dstrain_input,dstress_input,this->stress_shift,false,false);
//   EM sp0_dstrain = dstrain_input;

//   for (size_t i = 0 ; i < nstep ; i++ ) {
//     StateParameters sp(this->strain,this->stress,sp0_dstrain,dstress_input,this->stress_shift,false,false);

//     auto [p, R] = this->set_stress_variable(this->stress);
//     auto [dstrain,sp0] = this->plastic_deformation_strain(dstrain_input,dstress_input,sp);

//     this->stress += dstress_input;
//     this->strain += dstrain;
//     sp0_dstrain = sp0.dstrain;

//     auto [ev, gamma] = this->set_strain_variable(this->strain);
//     this->e = this->e0 - ev*(1.0+this->e0);
//   }
}

// ------------------------------------------------------------------- //
std::tuple<EM, EV, double> GHE_S::set_Dp_matrix(EV FEMdstrain) {
//   EM dstrain = this->FEMstrain_to_matrix(FEMdstrain);
//   EM dstress_input = EM::Zero(3,3);

//   StateParameters sp0(this->strain,this->stress,dstrain,dstress_input,this->stress_shift,false,false);
//   auto [ef1, ef2] = this->check_unload(sp0);

//   StateParameters sp(this->strain,this->stress,dstrain,dstress_input,this->stress_shift,ef1,ef2);
//   Eigen::Tensor<double,4> Ep = this->plastic_stiffness(sp);

//   EM dstress = this->solve_stress(dstrain,Ep);

//   StateParameters sp_check(this->strain,this->stress,dstrain,dstress,this->stress_shift,ef1,ef2);
//   auto [ef1_check, ef2_check] = this->check_unload(sp_check);

//   StateParameters sp2(strain,stress,dstrain,dstress,this->stress_shift,ef1_check,ef2_check);
//   this->update_parameters(sp2);

//   auto [ev, gamma] = this->set_strain_variable(this->strain);
//   this->e = this->e0 - ev*(1.0+this->e0);

//   this->stress += dstress;
//   this->strain += dstrain;

//   EM Dp = this->modulus_to_Dmatrix(Ep);
//   EV FEMstress = this->matrix_to_FEMstress(this->stress);
//   double stress_yy = -this->stress(1,1);
  Eigen::Tensor<double,4> Ep(3,3,3,3);
  Ep.setZero();
  EM Dp = this->modulus_to_Dmatrix(Ep);
  EV FEMstress = this->matrix_to_FEMstress(this->stress);
  double stress_yy = -this->stress(1,1);
  return {Dp, FEMstress, stress_yy};

}

// ------------------------------------------------------------------- //
EV GHE_S::strain_to_stress(EV FEMdstrain) {
//   EM dstrain = this->FEMstrain_to_matrix(FEMdstrain);
//   EM dstress_input = EM::Zero(3,3);

//   StateParameters sp(this->strain,this->stress,dstrain,dstress_input,this->stress_shift,false,false);
//   auto [p,R] = this->set_stress_variable(this->stress);

//   auto [dstress,_sp] = this->plastic_deformation_stress(dstrain,dstress_input,sp);

//   this->stress += dstress;
//   this->strain += dstrain;

//   auto [ev, gamma] = this->set_strain_variable(this->strain);
//   this->e = this->e0 - ev*(1.0+this->e0);

  EV FEMstress = this->matrix_to_FEMstress(this->stress);
  return FEMstress;
}
