#include "all.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include "elasto_plastic.h"
#include "ep_model.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;

// ------------------------------------------------------------------- //
StateParameters::StateParameters (EM strain, EM stress, EM dstrain, EM dstress, bool ef1, bool ef2)
  {
    this->strain = strain;
    this->stress = stress;
    this->dstrain = dstrain;
    this->dstress = dstress;
    this->pmin = 1.0;

    this->elastic_flag1 = ef1;
    this->elastic_flag2 = ef2;

    this->set_stress_variable();
    this->set_stress_increment();
  }

void StateParameters::set_stress_variable() {
  this->p = this->stress.trace()/3.0;
  this->sij = this->stress - this->p * EM::Identity(3,3);
  this->rij = this->sij / std::max(this->p,this->pmin);
  this->R = std::sqrt(1.5*this->rij.array().square().sum());
}

void StateParameters::set_stress_increment() {
  EM stress = this->stress + this->dstress;
  double p = stress.trace()/3.0;
  this->dp = p - this->p;
}

// ------------------------------------------------------------------- //
Li2002::Li2002 (double G0, double nu, double M, double eg, double d1,
                double c, double rlambdac, double xi,
                double m, double h1, double h2, double h3, double n,
                double d2, double h4, double a, double b1, double b2, double b3)
  {
    // Elastic parameters
    this->G0 = G0; this->nu = nu;
    // Critical state parameters
    this->M = M; this->c = c; this->eg = eg; this->rlambdac = rlambdac; this->xi = xi;
    // parameters associated with dr-mechanisms
    this->d1 = d1; this->m = m; this->h1 = h1; this->h2 = h2; this->h3 = h3; this->n = n;
    // parameters associated with dp-mechanisms
    this->d2 = d2; this->h4 = h4;
    // Default parameters
    this->a = a; this->b1 = b1; this->b2 = b2; this->b3 = b3;

    // minimum epsillon
    this->eps = 1.e-6;

    // stress parameters
    this->pr = 101.e3;
    this->pmin = 100.0;

    // stress and strain
    this->stress = EM::Zero(3,3);
    this->strain = EM::Zero(3,3);

    // BS parameters
    this->alpha = EM::Zero(3,3);
    this->beta = 0.0;
    this->H1 = 0.0;
    this->H2 = 0.0;

    // Accumulated index
    this->L1 = 0.0;

    // identity
    this->Z3 = EM::Zero(3,3);
    this->I3 = EM::Identity(3,3);

    // parameters
    this->sqrt2_3 = std::sqrt(2.0/3.0);
    this->fn = 2.0*(1+this->nu)/(3.0*(1.0-2.0*this->nu));
    this->rlambda_coeff = 2.0*this->nu/(1.0-2.0*this->nu);
    this->G2_coeff = this->fn*this->h4 / (this->h4 + this->sqrt2_3*this->fn*this->d2) / this->fn;
    this->g0 = this->c*(1.0+this->c) / (1.0+this->c*this->c);
    this->dg0 = (-this->c*this->c*(1-this->c)* std::pow(1+this->c,2)) / std::pow(1+this->c*this->c,3);
  }


// ------------------------------------------------------------------- //
void Li2002::print() {
  // std::cout << this->G0 << "\n";
}

void Li2002::clear_strain() {
  this->strain = EM::Zero(3,3);
}

// ------------------------------------------------------------------- //
void Li2002::initial_state(EV init_stress) {
  EM init_stress_mat = this->FEMstress_to_matrix(init_stress);

  // isotropic compression
  double compression_stress = init_stress_mat(0,0);
  this->isotropic_compression(this->e0,compression_stress);
  this->e = this->e0;

  // set initial parameters
  double p = this->set_stress_variable_p(this->stress);
  this->beta = p;
  this->H2 = p;

  size_t nstep = 50;
  EV dstrain_vec = EV::Zero(6);
  EV dstress_vec = EV::Zero(6);

  dstress_vec(2) = (init_stress_mat(2,2)-init_stress_mat(0,0))/nstep;

  EM dstrain_input = this->vector_to_matrix(dstrain_vec);
  EM dstress_input = this->vector_to_matrix(dstress_vec);

  StateParameters sp0(this->strain,this->stress,dstrain_input,dstress_input,false,false);

  for (size_t i = 0 ; i < nstep ; i++ ) {
    StateParameters sp(this->strain,this->stress,sp0.dstrain,dstress_input,false,false);

    auto [p, R] = this->set_stress_variable(this->stress);
    auto [dstrain,sp0] = this->plastic_deformation_strain(dstrain_input,dstress_input,sp);

    this->stress += dstress_input;
    this->strain += dstrain;

    auto [ev, gamma] = this->set_strain_variable(this->strain);
    this->e = this->e0 - ev*(1.0+this->e0);
  }
}

// ------------------------------------------------------------------- //
std::tuple<EM, EV> Li2002::set_Dp_matrix(EV FEMdstrain) {
  EM dstrain = this->FEMstrain_to_matrix(FEMdstrain);
  EM dstress_input = EM::Zero(3,3);

  StateParameters sp0(this->strain,this->stress,dstrain,dstress_input,false,false);
  auto [ef1, ef2] = this->check_unload(sp0);

  StateParameters sp(this->strain,this->stress,dstrain,dstress_input,ef1,ef2);
  Eigen::Tensor<double,4> Ep = this->plastic_stiffness(sp);

  EM dstress = this->solve_stress(dstrain,Ep);

  auto [ev, gamma] = this->set_strain_variable(this->strain);
  this->e = this->e0 - ev*(1.0+this->e0);

  this->stress += dstress;
  this->strain += dstrain;

  EM Dp = this->modulus_to_Dmatrix(Ep);
  EV FEMstress = this->matrix_to_FEMstress(this->stress);
  return {Dp, FEMstress};

}

// ------------------------------------------------------------------- //
EV Li2002::strain_to_stress(EV FEMdstrain) {
  EM dstrain = this->FEMstrain_to_matrix(FEMdstrain);
  EM dstress_input = EM::Zero(3,3);

  StateParameters sp(this->strain,this->stress,dstrain,dstress_input,false,false);
  auto [p,R] = this->set_stress_variable(this->stress);

  auto [dstress,_sp] = this->plastic_deformation_stress(dstrain,dstress_input,sp);

  this->stress += dstress;
  this->strain += dstrain;

  auto [ev, gamma] = this->set_strain_variable(this->strain);
  this->e = this->e0 - ev*(1.0+this->e0);

  EV FEMstress = this->matrix_to_FEMstress(this->stress);
  return FEMstress;
}

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Li2002::isotropic_compression(const double e0, const double compression_stress) {
  size_t nstep = 1000;
  double dcp = compression_stress / nstep;
  this->e = e0;

  EV dstress_vec = EV::Zero(6);
  dstress_vec(0) = dcp; dstress_vec(1) = dcp; dstress_vec(2) = dcp;
  EM dstress = EP::vector_to_matrix(dstress_vec);

  this->stress = EM::Zero(3,3);
  this->strain = EM::Zero(3,3);
  for (size_t i=0 ; i<nstep ; i++) {
    double p = this->set_stress_variable_p(this->stress);
    Eigen::Tensor<double,4> E = this->isotropic_compression_stiffness(this->e,p);
    EM dstrain = this->solve_strain(dstress,E);

    this->stress += dstress;
    this->strain += dstrain;
    double ev = this->set_strain_variable_ev(this->strain);
    this->e = e0 - ev*(1.0+e0);
  }

  this->clear_strain();
}

// ------------------------------------------------------------------- //
std::tuple<EM, StateParameters>
  Li2002::plastic_deformation_strain(const EM dstrain_given, const EM dstress_given, const StateParameters sp0) {
    auto [ef1, ef2] = this->check_unload(sp0);
    EM strain = sp0.strain;
    EM stress = sp0.stress;

    StateParameters sp(strain,stress,dstrain_given,dstress_given,ef1,ef2);
    Eigen::Tensor<double,4> Ep = this->plastic_stiffness(sp);
    EM dstrain_ep = this->solve_strain(dstress_given,Ep);

    StateParameters sp2(strain,stress,dstrain_ep,dstress_given,false,false);
    this->update_parameters(sp2);

    return {dstrain_ep, sp2};
}

// ------------------------------------------------------------------- //
std::tuple<EM, StateParameters>
  Li2002::plastic_deformation_stress(const EM dstrain_given, const EM dstress_given, const StateParameters sp0) {
    auto [ef1, ef2] = this->check_unload(sp0);
    EM strain = sp0.strain;
    EM stress = sp0.stress;

    StateParameters sp(strain,stress,dstrain_given,dstress_given,ef1,ef2);
    Eigen::Tensor<double,4> Ep = this->plastic_stiffness(sp);
    EM dstress_ep = this->solve_stress(dstrain_given,Ep);

    StateParameters sp2(strain,stress,dstrain_given,dstress_ep,false,false);
    this->update_parameters(sp2);

    return {dstress_ep, sp2};
}

// ------------------------------------------------------------------- //
std::tuple<double, double>
  Li2002::set_stress_variable(const EM stress) {
    double p = stress.trace()/3.0;
    EM r_stress = (stress - p*EM::Identity(3,3)) / std::max(p,this->pmin);
    double R = std::sqrt(1.5*r_stress.array().square().sum());
    return {p, R};
  }

double Li2002::set_stress_variable_p(const EM stress) {
  double p = stress.trace()/3.0;
  return p;
}

std::tuple<double,double>
  Li2002::set_strain_variable(const EM strain) {
    double ev = strain.trace();
    EM dev_strain = strain - ev/3.0 * this->I3;
    double gamma = std::sqrt(2.0/3.0)*dev_strain.norm();
    return {ev, gamma};
  }

double Li2002::set_strain_variable_ev(const EM strain) {
  double ev = strain.trace();
  return ev;
}

// ------------------------------------------------------------------- //
Eigen::Tensor<double,4>
  Li2002::isotropic_compression_stiffness(const double e, const double p) {
    double G = this->elastic_modulus_G(e,p);
    double G2 = G*this->G2_coeff;
    Eigen::Tensor<double,4> E2 = this->elastic_stiffness(G2);
    return E2;
  }

Eigen::Tensor<double,4> Li2002::elastic_stiffness(const double G) {
  double rlambda = G*this->rlambda_coeff;
  Eigen::Tensor<double,4> Dijkl = this->set_Dijkl();
  Eigen::Tensor<double,4> Dikjl = this->set_Dikjl();
  Eigen::Tensor<double,4> Ee = rlambda*Dijkl + 2*G*Dikjl;
  return Ee;
}

double Li2002::elastic_modulus_G(const double e, const double p) {
  double G = this->G0*std::pow(2.97-e,2) / (1+e) * std::sqrt(std::max(p,this->pmin)*this->pr);
  return G;
}

std::tuple<double, double>
  Li2002::elastic_modulus(const double e, const double p) {
    double G = this->G0*std::pow(2.97-e,2) / (1+e) * std::sqrt(std::max(p,this->pmin)*this->pr);
    double K = G*2*(1+this->nu)/(3*(1-2*this->nu));
    return {G, K};
  }

std::tuple<double, double>
  Li2002::elastic_modulus_lame() {
    double p = this->set_stress_variable_p(this->stress);
    auto [G0, K0] = this->elastic_modulus(this->e,p);
    return {G0, K0 - G0*2.0/3.0};
  }

// ------------------------------------------------------------------- //
EM Li2002::solve_strain(const EM stress_mat, const Eigen::Tensor<double,4> E) {
  EV b = EV::Zero(9);
  b(0) = stress_mat(0,0); b(1) = stress_mat(0,1); b(2) = stress_mat(0,2);
  b(3) = stress_mat(1,0); b(4) = stress_mat(1,1); b(5) = stress_mat(1,2);
  b(6) = stress_mat(2,0); b(7) = stress_mat(2,1); b(8) = stress_mat(2,2);

  EM A = EM::Zero(9,9);
  for(int i=0;i<3;i++){ for(int j=0;j<3;j++) {
    size_t i0 = 3*i+j;
    for(int k=0;k<3;k++){ for(int l=0;l<3;l++) {
      size_t j0 = 3*k+l;
      A(i0,j0) = E(i,j,k,l);
    }}
  }}

  EV x = A.partialPivLu().solve(b);

  EM strain_mat = EM::Zero(3,3);
  strain_mat(0,0) = x(0); strain_mat(0,1) = x(1); strain_mat(0,2) = x(2);
  strain_mat(1,0) = x(3); strain_mat(1,1) = x(4); strain_mat(1,2) = x(5);
  strain_mat(2,0) = x(6); strain_mat(2,1) = x(7); strain_mat(2,2) = x(8);

  return strain_mat;
}

// ------------------------------------------------------------------- //
EM Li2002::solve_stress(const EM strain_mat, const Eigen::Tensor<double,4> E) {
  EM stress_mat = EM::Zero(3,3);
  for(int i=0;i<3;i++){ for(int j=0;j<3;j++) {
    for(int k=0;k<3;k++){ for(int l=0;l<3;l++) {
      stress_mat(i,j) = stress_mat(i,j) + E(i,j,k,l) * strain_mat(k,l);
    }}
  }}

  return stress_mat;
}

// ------------------------------------------------------------------- //
std::tuple<double, double>
  Li2002::check_unload(StateParameters sp) {
    this->set_mapping_stress(sp);
    this->set_parameters(sp);
    this->set_parameter_nm(sp);
    this->set_parameter_TZ(sp);

    bool elastic_flag1,elastic_flag2;

    double dL1 = (sp.Tij.array() * sp.dstrain.array()).sum();
    if (dL1 < 0.0) {
      elastic_flag1 = true;
      this->alpha = sp.rij;
    } else {
      elastic_flag1 = false;
    }

    double dL2 = (sp.Zij.array() * sp.dstrain.array()).sum();
    if (dL2 < 0.0) {
      elastic_flag2 = true;
      this->beta = sp.p;
    } else {
      elastic_flag2 = false;
    }

    return {elastic_flag1, elastic_flag2};
  }

// ------------------------------------------------------------------- //
Eigen::Tensor<double,4>
  Li2002::plastic_stiffness(StateParameters &sp) {
    this->set_mapping_stress(sp);
    this->set_parameters(sp);
    this->set_parameter_nm(sp);
    this->set_parameter_TZ(sp);
    Eigen::Tensor<double,4> Ep = this->set_tensor_Ep(sp);
    return Ep;
  }

// ------------------------------------------------------------------- //
void Li2002::update_parameters(StateParameters &sp) {
  this->set_mapping_stress(sp);
  this->set_parameters(sp);
  this->set_parameter_nm(sp);
  this->set_parameter_TZ(sp);
  double dL1 = (sp.Tij.array() * sp.dstrain.array()).sum();
  double dL2 = (sp.Zij.array() * sp.dstrain.array()).sum();
  if (!sp.elastic_flag1) {
    this->L1 += dL1;
    this->H1 += sp.Kp1_b*dL1 / sp.p;
  }
  if (!sp.elastic_flag2) {
    this->H2 += sp.Kp2_b*dL2;
  }
}

// ------------------------------------------------------------------- //
void Li2002::set_mapping_stress(StateParameters &sp) {

    if (sp.elastic_flag1) {
      this->alpha = sp.rij;
    }
    if (sp.elastic_flag2) {
      this->beta = sp.p;
    }

    if ((sp.rij-this->alpha).norm() < 1.e-6) {
      sp.elastic_flag1 = true;
    } else {
      auto [F1,rij_bar,R_bar,g_bar] = this->_F1_boundary_surface_all(1.0,sp.rij,this->alpha);
      if (F1 > 0.0) {
        sp.rij_bar = rij_bar; sp.R_bar = R_bar; sp.g_bar = g_bar;
        this->H1 = sp.R_bar/sp.g_bar;
        sp.rho1_ratio = 1.0;
      } else {
        double t = this->_find_rho1_ratio(sp.rij,this->alpha);
        auto [rij_bar, R_bar, g_bar] = this->_mapping_r(t,sp.rij,this->alpha);
        sp.rij_bar = rij_bar; sp.R_bar = R_bar; sp.g_bar = g_bar;
        sp.rho1_ratio = t;
      }
    }

    if (std::abs(sp.p-this->beta) == 0.0) {
      sp.elastic_flag2 = true;
    } else {
      if (sp.p > this->H2) {
        this->H2 = sp.p;
      }
      if (sp.dp > 0.0) {
        if (sp.p <= this->beta) {
          sp.elastic_flag2 = true;
          return;
        }
        sp.p_bar = this->H2;
      } else if (sp.dp < 0.0) {
        if (this->beta <= sp.p) {
          sp.elastic_flag2 = true;
          return;
        }
        sp.p_bar = this->pmin;
      } else {
        sp.elastic_flag2 = true;
        return;
      }
      double rho2 = std::abs(sp.p - this->beta);
      double rho2_b = std::abs(sp.p_bar - this->beta);
      sp.rho2_ratio = rho2_b / rho2;
    }
}

std::tuple<double, EM, double, double>
  Li2002::_F1_boundary_surface_all(const double t, const EM rij, const EM alpha) {
    auto [rij_bar, R_bar, g_bar] = this->_mapping_r(t,rij,alpha);
    return {R_bar-this->H1*g_bar, rij_bar, R_bar, g_bar};
  }

double Li2002::_F1_boundary_surface(const double t, const EM rij, const EM alpha) {
  auto [R_bar, g_bar] = this->_mapping_r_Rg(t,rij,alpha);
  return R_bar-this->H1*g_bar;
}

std::tuple<EM, double, double>
  Li2002::_mapping_r(const double t, const EM rij, const EM alpha) {
    EM rij_bar = alpha + t*(rij-alpha);
    auto [g_bar, J2] = this->g_theta_J2(rij_bar);
    double R_bar = std::sqrt(3.0*J2);
    return {rij_bar, R_bar, g_bar};
  }

std::tuple<double, double>
  Li2002::_mapping_r_Rg(const double t, const EM rij, const EM alpha) {
    EM rij_bar = alpha + t*(rij-alpha);
    auto [g_bar, J2] = this->g_theta_J2(rij_bar);
    double R_bar = std::sqrt(3.0*J2);
    return {R_bar, g_bar};
  }

double Li2002::_find_rho1_ratio(const EM rij, const EM alpha) {
  int nstep = 100;
  double xtol = 1.e-8;
  double rtol = 1.e-10;
  double gr = (std::sqrt(5.0) + 1.0)/2.0;

  double t0 = 0.99;
  double t1 = 1.e2;

  double F1_t0 = this->_F1_boundary_surface(t0,rij,alpha);
  double F1_t1 = this->_F1_boundary_surface(t1,rij,alpha);

  for (size_t i = 0 ; i < nstep ; i++ ){
    if (F1_t0 * F1_t1 >= 0.0) {
      // t0 = 0.1*t0;
      t1 = 10.0*t1;
      // F1_t0 = this->_F1_boundary_surface(t0,rij,alpha);
      F1_t1 = this->_F1_boundary_surface(t1,rij,alpha);
      // std::cout << t0 << ": "<< F1_t0 << std::endl;
      // std::cout << t1 << ": "<< F1_t1 << std::endl;
    } else {
      break;
    }
  }

  double tl = t1 - (t1-t0)/gr;
  double tr = t0 + (t1-t0)/gr;
  double F1_tl = this->_F1_boundary_surface(tl,rij,alpha);
  double F1_tr = this->_F1_boundary_surface(tr,rij,alpha);

  for (size_t i = 0 ; i < nstep ; i++ ){
    if (std::abs(tl-tr) < xtol) {
      break;
    }
    if (std::abs(F1_tl-F1_tr) < rtol) {
      break;
    }
    if (std::abs(F1_tl) < std::abs(F1_tr)) {
      t1 = tr;
    } else {
      t0 = tl;
    }
    tl = t1 - (t1-t0)/gr;
    tr = t0 + (t1-t0)/gr;
    F1_tl = this->_F1_boundary_surface(tl,rij,alpha);
    F1_tr = this->_F1_boundary_surface(tr,rij,alpha);
    // std::cout << i << std::endl;
    // std::cout << tl << ": "<< F1_tl << std::endl;
    // std::cout << tr << ": "<< F1_tr << std::endl;
  }
  return (tl+tr)/2.0;
}

// ------------------------------------------------------------------- //
void Li2002::set_parameters(StateParameters &sp) {
  auto [Ge,Ke] = this->elastic_modulus(this->e,sp.p);
  sp.Ge = Ge; sp.Ke = Ke;
  sp.psi = this->state_parameter(this->e,sp.p);
  sp.g = this->g_theta(sp.sij);

  if (sp.elastic_flag1) {
    sp.Kp1_b = 0.0;
    sp.D1 = 0.0;
  } else {
    double h = this->_scaling_factor(this->e,sp.rho1_ratio);
    auto [Kp1,Kp1_b] = this->_plastic_modulus1(sp.Ge,sp.R_bar,sp.g_bar,sp.rho1_ratio,h,sp.psi);
    sp.h = h; sp.Kp1 = Kp1; sp.Kp1_b = Kp1_b;
    sp.D1 = this->_dilatancy1(sp.R,sp.g,sp.rho1_ratio,sp.psi);
  }

  if (sp.elastic_flag2 || (sp.R == 0.0)) {
    sp.Kp2_b = 0.0;
    sp.D2 = 0.0;
  } else {
    double sign = sp.dp / std::abs(sp.dp);
    double Mg_R = this->M*sp.g/sp.R;
    auto [Kp2,Kp2_b] = this->_plastic_modulus2(sp.Ge,Mg_R,sp.rho2_ratio,sign);
    sp.Kp2 = Kp2; sp.Kp2_b = Kp2_b;
    sp.D2 = this->_dilatancy2(Mg_R,sign);
  }
}

double Li2002::_accumulated_load_index(const double L1) {
  double fL = 1.0-this->b3;
  double fL2 = std::pow(1.0-L1/this->b1,2)+(L1/this->b1)/std::pow(this->b2,2) + this->b3;
  return fL/std::sqrt(fL2);
}

double Li2002::_scaling_factor(const double e, const double rho1_ratio) {
  double fL = this->_accumulated_load_index(this->L1);
  double r1 = std::pow(1.0/rho1_ratio,10);
  double h = (this->h1-this->h2*e)*(r1+this->h3*fL*(1.0-r1));
  return h;
}

std::tuple<double, double>
  Li2002::_plastic_modulus1(const double G, const double R_bar, const double g_bar, const double rho1_ratio, const double h, const double psi) {
    double Mg_R = this->M*g_bar*std::exp(-this->n*psi) / R_bar;
    double Kp1 = G*h*(Mg_R*rho1_ratio - 1.0);
    double Kp1_b = G*h*(Mg_R - 1.0);
    return {Kp1, Kp1_b};
  }

std::tuple<double, double>
  Li2002::_plastic_modulus2(const double G, const double Mg_R, const double rho2_ratio, const double sign) {
    double Kp2 = G*this->h4*Mg_R*sign * std::pow(rho2_ratio,this->a);
    if ((rho2_ratio == 1.0) && sign > 0.0) {
      return {Kp2, Kp2};
    } else {
      return {Kp2, 0.0};
    }
  }

double Li2002::_dilatancy1(const double R, const double g, const double rho1_ratio, const double psi) {
  double R_Mg = R / (this->M*g);
  double D1 = this->d1*(std::exp(this->m*psi)*std::sqrt(rho1_ratio) - R_Mg);
  return D1;
}

double Li2002::_dilatancy2(const double Mg_R, const double sign) {
  if (Mg_R >= 1.0) {
    return this->d2*(Mg_R-1.0)*sign;
  } else {
    return 0.0;
  }
}

// ------------------------------------------------------------------- //
void Li2002::set_parameter_nm(StateParameters &sp) {
  if (!sp.elastic_flag1) {
    double theta3_bar = this->Lode_angle(sp.sij);
    double dg_bar = this->dg_theta(theta3_bar);
    EM dF1 = this->_dF1_r(sp.rij_bar,sp.R_bar,theta3_bar,sp.g_bar,dg_bar);
    double dF1_tr = dF1.trace();
    EM nij = dF1 - this->I3 * dF1_tr/3.0;
    sp.nij = nij / nij.norm();
  }

  double r_abs = sp.rij.norm();
  if (r_abs == 0.0) {
    sp.mij = this->Z3;
  } else {
    sp.mij = sp.rij / r_abs;
  }
}

EM Li2002::_dF1_r(const EM rij_bar, const double R_bar, const double theta3_bar, const double g_bar, const double dg_bar) {
  if (std::abs(R_bar) < this->eps) {
    return this->Z3;
  }
  double st_bar = std::sin(theta3_bar);
  double a = R_bar*g_bar + 3*R_bar*st_bar*dg_bar;
  double b = 9*dg_bar;
  double c = 1.5/std::pow(R_bar*g_bar,2);
  EM rr_bar = rij_bar * rij_bar.transpose();
  return (a*rij_bar + b*rr_bar)*c;
}

// ------------------------------------------------------------------- //
void Li2002::set_parameter_TZ(StateParameters &sp) {
  if (sp.elastic_flag1 || sp.elastic_flag2 || (sp.R == 0.0)) {
    sp.B = 0.0;
  } else {
    double nm = (sp.nij.array() * sp.mij.array()).sum();
    double nr = (sp.nij.array() * sp.rij.array()).sum();
    double Bu = 2*sp.Ge*nm - this->sqrt2_3*sp.Ke*sp.D2*nr;
    double Bd = this->sqrt2_3*sp.Ke*sp.D2 + sp.Kp2;
    sp.B = Bu / Bd;
  }

  if (sp.elastic_flag1) {
    sp.Tij = this->Z3;
  } else {
    double nr = (sp.nij.array() * sp.rij.array()).sum();
    EM Tu = 2*sp.Ge*sp.nij - sp.Ke*(nr+sp.B)*this->I3;
    double Td = 2*sp.Ge - this->sqrt2_3*sp.Ke*sp.D1*(nr+sp.B) + sp.Kp1;
    sp.Tij = Tu / Td;
  }

  if (sp.elastic_flag2) {
    sp.Zij = this->Z3;
  } else {
    EM Zu = sp.Ke*this->I3 - this->sqrt2_3*sp.Ke*sp.D1*sp.Tij;
    double Zd;
    if (sp.R == 0.0) {
      double Kp2_D2 = sp.Ge*this->h4/this->d2 * sp.rho2_ratio*this->a;
      Zd = this->sqrt2_3*sp.Ke + Kp2_D2;
    } else {
      Zd = this->sqrt2_3*sp.Ke*sp.D2 + sp.Kp2;
    }
    sp.Zij = Zu / Zd;
  }

}

// ------------------------------------------------------------------- //
Eigen::Tensor<double,4>
  Li2002::set_tensor_Ep(StateParameters &sp) {
    Eigen::Tensor<double,4> Lm0 = this->set_Dikjl();
    Eigen::Tensor<double,4> Lm1(3,3,3,3), Lm2(3,3,3,3);

    if (sp.elastic_flag1) {
      Lm1.setZero();
    } else {
      EM nD = sp.nij + std::sqrt(2.0/27.0)*sp.D1*this->I3;
      for(int i=0;i<3;i++){ for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++){ for(int l=0;l<3;l++){
          Lm1(i,j,k,l) = nD(i,j) * sp.Tij(k,l);
        }}
      }}
    }

    if (sp.elastic_flag2) {
      Lm2.setZero();
    } else if (sp.R == 0.0){
      EM mD = std::sqrt(2.0/27.0)*this->I3;
      for(int i=0;i<3;i++){ for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++){ for(int l=0;l<3;l++){
          Lm2(i,j,k,l) = mD(i,j) * sp.Zij(k,l);
        }}
      }}
    } else {
      EM mD = sp.mij + std::sqrt(2.0/27.0)*sp.D2*this->I3;
      for(int i=0;i<3;i++){ for(int j=0;j<3;j++) {
        for(int k=0;k<3;k++){ for(int l=0;l<3;l++){
          Lm2(i,j,k,l) = mD(i,j) * sp.Zij(k,l);
        }}
      }}
    }

    Eigen::Tensor<double,4> Lm = Lm0 - Lm1 - Lm2;
    Eigen::Tensor<double,4> Ee = this->elastic_stiffness(sp.Ge);

    Eigen::Tensor<double,4> Ep(3,3,3,3);
    Ep.setZero();

    for(int i=0;i<3;i++){ for(int j=0;j<3;j++) {
      for(int k=0;k<3;k++){ for(int l=0;l<3;l++){
        for(int p=0;p<3;p++){ for(int q=0;q<3;q++){
          Ep(i,j,k,l) = Ep(i,j,k,l) + Ee(i,j,p,q)*Lm(p,q,k,l);
        }}
      }}
    }}

    return Ep;
  }

// ------------------------------------------------------------------- //
double Li2002::state_parameter(const double e, const double p) {
  double psi = e - (this->eg-this->rlambdac* std::pow(std::max(p,this->pmin)/this->pr,this->xi));
  return psi;
}

double Li2002::Lode_angle(const EM dev_stress) {
  double J2 = 0.5*dev_stress.array().square().sum();
  if (J2 == 0.0) {
    return 0.0;
  }
  double J3 = -dev_stress.determinant();
  double s3 = 0.5*J3 * std::pow(3.0/J2,1.5);
  s3 = std::max(s3,-1.0);
  s3 = std::min(s3,1.0);
  double theta3 = std::asin(s3);
  return theta3;
}

std::tuple<double, double>
  Li2002::Lode_angle_J2(const EM dev_stress) {
    double J2 = 0.5*dev_stress.array().square().sum();
    if (J2 == 0.0) {
      return {0.0, 0.0};
    }
    double J3 = -dev_stress.determinant();
    double s3 = 0.5*J3 * std::pow(3.0/J2,1.5);
    s3 = std::max(s3,-1.0);
    s3 = std::min(s3,1.0);
    double theta3 = std::asin(s3);
    return {theta3, J2};
  }

double Li2002::g_theta(const EM dev_stress) {
  double theta3 = this->Lode_angle(dev_stress);
  if (theta3 == 0.0) {
    return 0.0;
  }
  double st = std::sin(theta3);
  double g1 = std::sqrt(std::pow(1.0+this->c*this->c,2) + 4*this->c*(1.0-this->c*this->c)*st) - (1.0+this->c*this->c);
  double g2 = 2.0*(1.0-this->c)*st;
  return g1/g2;
}

std::tuple<double, double>
  Li2002::g_theta_J2(const EM dev_stress) {
    auto [theta3,J2] = this->Lode_angle_J2(dev_stress);
    if (theta3 == 0.0) {
      return {this->g0, J2};
    }
    double st = std::sin(theta3);
    double g1 = std::sqrt(std::pow(1.0+this->c*this->c,2) + 4*this->c*(1.0-this->c*this->c)*st) - (1.0+this->c*this->c);
    double g2 = 2.0*(1.0-this->c)*st;
    return {g1/g2, J2};
  }

double Li2002::dg_theta(const double theta3) {
  double st = std::sin(theta3);
  if (st == 0.0) {
    return this->dg0;
  }
  double g11 = this->c*(1+this->c);
  double g12 = st*std::sqrt(std::pow(1+this->c*this->c,2) + 4*this->c*(1-this->c*this->c)*st);
  double g21 = std::sqrt(std::pow(1+this->c*this->c,2)+4*this->c*(1-this->c*this->c)*st) - (1+this->c*this->c);
  double g22 = 2*(1-this->c)*st*st;
  double dg = g11/g12 - g21/g22;
  return dg;
}

// ------------------------------------------------------------------- //
Eigen::Tensor<double,4> Li2002::set_Dijkl() {
    Eigen::Tensor<double,4> Dijkl(3,3,3,3);
    Dijkl.setZero();

    for (int i=0;i<3;i++){
      for (int k=0;k<3;k++){
        Dijkl(i,i,k,k) = 1.0;
      }
    }
    return Dijkl;
  }

Eigen::Tensor<double,4> Li2002::set_Dikjl() {
    Eigen::Tensor<double,4> Dikjl(3,3,3,3);
    Dikjl.setZero();

    for (int i=0;i<3;i++){
      for (int k=0;k<3;k++){
        Dikjl(i,k,i,k) = 1.0;
      }
    }
    return Dikjl;
  }
