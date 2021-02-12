#include "all.h"
#include <Eigen/Dense>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"


// ------------------------------------------------------------------- //
Element::Element (size_t id, std::string style, int material_id, std::vector<size_t> inode)
  {
    Element::id = id;
    Element::style = style;
    Element::material_id = material_id;
    Element::inode = inode;

    Element::gravity = 9.8;
    Element::nnode = inode.size();

    Element::set_style();
  }

// ------------------------------------------------------------------- //
void
  Element::print() {
    std::cout << Element::id << ": ";
    std::cout << Element::style << ", ";
    std::cout << Element::material_id << ", ";

    for (size_t i = 0 ; i < Element::inode.size() ; i++) {
      std::cout << Element::inode.at(i) << " ";
    }
    std::cout << "\n";
  }

void
  Element::set_style() {
    ElementStyle* estyle_p = set_element_style(Element::style);

    Element::dim = estyle_p->dim;
    Element::ng = estyle_p->ng;
    Element::xi = estyle_p->xi;
    Element::w  = estyle_p->w;

    Element::ng_all = estyle_p->ng_all;
    Element::n_list  = estyle_p->n_list;
    Element::dn_list = estyle_p->dn_list;
    Element::w_list = estyle_p->w_list;

    Element::dn_center = estyle_p->dn_center;
  }

void
  Element::set_nodes(std::vector<Node*> nodes_p) {
    Element::nodes_p = nodes_p;
  }

void
  Element::set_material(Material* material_p) {
    if (material_p == nullptr) {
      Element::material_p = nullptr;
      Element::rho = 0.0;
    } else {
      Element::material_p = material_p;
      Element::rho = material_p->rho;
    }
  }

void
  Element::set_pointer_list(){
    Element::u_p.resize(Element::nnode);
    Element::v_p.resize(Element::nnode);

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      Element::u_p[inode] = &Element::nodes_p[inode]->u;
      Element::v_p[inode] = &Element::nodes_p[inode]->v;
    }
  }

void
  Element::set_xn(){
    Element::xnT = Eigen::MatrixXd::Zero(2,Element::nnode);

    for (size_t inode = 0 ; inode < Element::nnode ; inode++ ) {
      Element::xnT(0,inode) = Element::nodes_p[inode]->xyz[0] + (*Element::u_p[inode])[0];
      Element::xnT(1,inode) = Element::nodes_p[inode]->xyz[1] + (*Element::u_p[inode])[1];
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_local_matrix_init(const size_t dof){
    Element::dof = dof;
    Element::ndof = dof*Element::nnode;

    Element::M_diag = Eigen::VectorXd::Zero(Element::ndof);

    Element::K = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);
    Element::K_diag = Eigen::VectorXd::Zero(Element::ndof);
    Element::K_off_diag = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);

    Element::C = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);
    Element::C_diag = Eigen::VectorXd::Zero(Element::ndof);
    Element::C_off_diag = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);

    Element::force = Eigen::VectorXd::Zero(Element::ndof);

    if (Element::dim == 2) {
      double V = 0.0;
      for (size_t i = 0 ; i < Element::ng_all ; i++){
        double det;
        Eigen::Matrix2d jacobi;

        std::tie(det, jacobi) = Element::mk_jacobi(Element::xnT, Element::dn_list[i]);

        V += det * Element::w_list[i];
      }
      Element::mass = Element::rho * V;

    } else if (Element::dim == 1) {
      Element::imp = Element::material_p->mk_imp(Element::dof);
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_local_matrix() {
    if (Element::dim == 2) {
      Eigen::MatrixXd M = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);

      Element::C = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);
      Element::K = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);

      Element::De = Element::material_p->mk_d(Element::dof);
      Element::Dv = Element::material_p->mk_visco(Element::dof);

      for (size_t i = 0 ; i < Element::ng_all ; i++){
        double det, detJ;
        Eigen::MatrixXd dnj, N, Me, B, K, C;

        std::tie(det, dnj) = Element::mk_dnj(Element::xnT, Element::dn_list[i]);

        N = Element::mk_n(Element::dof, Element::nnode, Element::n_list[i]);
        Me = mk_m(N);

        B = Element::mk_b(Element::dof, Element::nnode, dnj);
        K = Element::mk_k(B, Element::De);
        C = Element::mk_k(B, Element::Dv);

        detJ = det * Element::w_list[i];

        M += Me * detJ;
        Element::C += C * detJ;
        Element::K += K * detJ;
      }

      double tr_M = M.trace() / Element::dof;
      Element::M_diag = M.diagonal() * Element::mass/tr_M;

      Element::K_diag = Element::K.diagonal();
      Element::K_off_diag = Element::K_diag.asDiagonal();
      Element::K_off_diag = Element::K - Element::K_off_diag;

      Element::C_diag = Element::C.diagonal();
      Element::C_off_diag = Element::C_diag.asDiagonal();
      Element::C_off_diag = Element::C - Element::C_off_diag;

    } else if (Element::dim == 1) {
      if (Element::style.find("input") != std::string::npos) {
        Element::C = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);

        for (size_t i = 0 ; i < Element::ng_all ; i++){
          double det, detJ;
          Eigen::MatrixXd N, q, NqN;

          std::tie(det, q) =
              Element::mk_q(Element::dof, Element::xnT,  Element::dn_list[i]);

          N = Element::mk_n(Element::dof, Element::nnode, Element::n_list[i]);
          NqN = Element::mk_nqn(N, q, Element::imp);

          detJ = det * Element::w_list[i];
          Element::C += NqN * detJ;
        }

        Element::C_diag = Element::C.diagonal();
        Element::C_off_diag = Element::C_diag.asDiagonal();
        Element::C_off_diag = Element::C - Element::C_off_diag;
      }
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_local_vector() {
    if (Element::dim == 2) {
      Element::force = Eigen::VectorXd::Zero(Element::ndof);

      double V = 0.0;
      for (size_t i = 0 ; i < Element::ng_all ; i++){
        double det, detJ;
        Eigen::Matrix2d jacobi;
        Eigen::MatrixXd N;

        std::tie(det, jacobi) = Element::mk_jacobi(Element::xnT, Element::dn_list[i]);
        N = Element::mk_n(Element::dof, Element::nnode, Element::n_list[i]);

        detJ = det * Element::w_list[i];

        V += detJ;
        Element::force += N.row(1)*detJ * Element::gravity;
      }

      Element::force *= Element::mass / V;
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_local_update() {
    if (Element::dim == 2) {
      Eigen::MatrixXd M = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);
      Element::C = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);
      Element::Dv = Element::material_p->mk_visco(Element::dof);
      Element::force = Eigen::VectorXd::Zero(Element::ndof);

      double V = 0.0;
      for (size_t i = 0 ; i < Element::ng_all ; i++){
        double det, detJ;
        Eigen::MatrixXd dnj, N, Me, B, K, C;

        std::tie(det, dnj) = Element::mk_dnj(Element::xnT, Element::dn_list[i]);

        N = Element::mk_n(Element::dof, Element::nnode, Element::n_list[i]);
        Me = mk_m(N);

        B = Element::mk_b(Element::dof, Element::nnode, dnj);
        C = Element::mk_k(B, Element::Dv);

        detJ = det * Element::w_list[i];
        M += Me * detJ;
        Element::C += C * detJ;

        V += detJ;
        Element::force += N.row(1)*detJ * Element::gravity;
      }

      double tr_M = M.trace() / Element::dof;
      Element::M_diag = M.diagonal() * Element::mass/tr_M;

      Element::C_diag = Element::C.diagonal();
      Element::C_off_diag = Element::C_diag.asDiagonal();
      Element::C_off_diag = Element::C - Element::C_off_diag;

      Element::force *= Element::mass / V;

    } else if (Element::dim == 1) {
      if (Element::style.find("input") != std::string::npos) {
        Element::C = Eigen::MatrixXd::Zero(Element::ndof,Element::ndof);

        for (size_t i = 0 ; i < Element::ng_all ; i++){
          double det, detJ;
          Eigen::MatrixXd N, q, NqN;

          std::tie(det, q) =
              Element::mk_q(Element::dof, Element::xnT,  Element::dn_list[i]);

          N = Element::mk_n(Element::dof, Element::nnode, Element::n_list[i]);
          NqN = Element::mk_nqn(N, q, Element::imp);

          detJ = det * Element::w_list[i];
          Element::C += NqN * detJ;
        }

        Element::C_diag = Element::C.diagonal();
        Element::C_off_diag = Element::C_diag.asDiagonal();
        Element::C_off_diag = Element::C - Element::C_off_diag;
      }
    }
  }


// ------------------------------------------------------------------- //
void
  Element::mk_ku() {
    Eigen::VectorXd u, ku;

    u = Element::mk_u_hstack(Element::u_p);
    ku = Element::K * u;

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      size_t i0 = inode*Element::dof;
      for (size_t i = 0 ; i < Element::dof ; i++) {
        Element::nodes_p[inode]->force(i) += ku(i0+i);
      }
    }
  }

void
  Element::mk_ku_u(const std::vector<Eigen::VectorXd*> u_p) {
    Eigen::VectorXd u, ku;

    u = Element::mk_u_hstack(u_p);
    ku = Element::K * u;

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      size_t i0 = inode*Element::dof;
      for (size_t i = 0 ; i < Element::dof ; i++) {
        Element::nodes_p[inode]->force(i) += ku(i0+i);
      }
    }

  }

void
  Element::mk_cv() {
    Eigen::VectorXd v, cv;

    v = Element::mk_u_hstack(Element::v_p);
    cv = Element::C_off_diag * v;

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      size_t i0 = inode*Element::dof;
      for (size_t i = 0 ; i < Element::dof ; i++) {
        Element::nodes_p[inode]->force(i) += cv(i0+i);
      }
    }
  }

void
  Element::mk_ku_cv() {
    Eigen::VectorXd u, ku;
    Eigen::VectorXd v, cv;

    u = Element::mk_u_hstack(Element::u_p);
    v = Element::mk_u_hstack(Element::v_p);
    ku = Element::K * u;
    cv = Element::C_off_diag * v;

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      size_t i0 = inode*Element::dof;
      for (size_t i = 0 ; i < Element::dof ; i++) {
        Element::nodes_p[inode]->force(i) += ku(i0+i) + cv(i0+i);
      }
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_B_stress() {
    if (Element::dim == 1) {
      Element::mk_ku();

    } else if (Element::dim == 2) {
      Eigen::VectorXd force = Eigen::VectorXd::Zero(Element::ndof);
      Eigen::MatrixXd u = Element::mk_u_vstack(Element::u_p);

      for (size_t i = 0 ; i < Element::ng_all ; i++){
        double det, detJ;
        Eigen::MatrixXd dnj, BT;
        Eigen::VectorXd stress;

        std::tie(det, dnj) = Element::mk_dnj(Element::xnT, Element::dn_list[i]);
        BT = Element::mk_b_T(Element::dof, Element::nnode, dnj);
        stress = Element::Hencky_stress(Element::De, dnj, u);

        detJ = det * Element::w_list[i];
        force += BT * stress * detJ;
      }

      for (size_t inode = 0 ; inode < Element::nnode ; inode++){
        size_t i0 = inode*Element::dof;
        for (size_t i = 0 ; i < Element::dof ; i++) {
          Element::nodes_p[inode]->force(i) += force(i0+i);
        }
      }
    }

  }


// ------------------------------------------------------------------- //
Eigen::VectorXd
  Element::mk_u_hstack(const std::vector<Eigen::VectorXd*> u_p) {
    Eigen::VectorXd u(Element::ndof);

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      size_t i0 = inode*Element::dof;
      for (size_t i = 0 ; i < Element::dof ; i++) {
        u(i0+i) = (*u_p[inode])(i);
      }
    }
    return u;
  }

Eigen::MatrixXd
  Element::mk_u_vstack(const std::vector<Eigen::VectorXd*> u_p) {
    Eigen::MatrixXd u(Element::nnode,Element::dof);

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      for (size_t i = 0 ; i < Element::dof ; i++) {
        u(inode,i) = (*u_p[inode])(i);
      }
    }
    return u;
  }


// ------------------------------------------------------------------- //
void
  Element::update_inputwave(const Eigen::VectorXd vel0) {
    Eigen::VectorXd v(Element::ndof), cv(Element::ndof);
    // Eigen::VectorXd cv(Element::ndof);

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      size_t i0 = inode*Element::dof;
      for (size_t i = 0 ; i < Element::dof ; i++) {
        v(i0+i) = vel0(i);
      }
    }
    cv = Element::C * v;

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      size_t i0 = inode*Element::dof;
      for (size_t i = 0 ; i < Element::dof ; i++) {
        Element::nodes_p[inode]->force(i) -= 2.0*cv(i0+i);
      }
    }
  }

// ------------------------------------------------------------------- //
void
  Element::update_bodyforce(const Eigen::VectorXd acc0) {
    Element::mk_bodyforce(acc0);

    for (size_t inode = 0 ; inode < Element::nnode ; inode++){
      size_t i0 = inode*Element::dof;
      for (size_t i = 0 ; i < Element::dof ; i++) {
        Element::nodes_p[inode]->force(i) -= Element::force(i0+i);
      }
    }
  }

void
  Element::mk_bodyforce(const Eigen::VectorXd acc0) {
    if (Element::dim == 2) {
      Element::force = Eigen::VectorXd::Zero(Element::ndof);

      double V = 0.0;
      for (size_t i = 0 ; i < Element::ng_all ; i++){
        double det, detJ;
        Eigen::Matrix2d jacobi;
        Eigen::MatrixXd N;

        std::tie(det, jacobi) = Element::mk_jacobi(Element::xnT, Element::dn_list[i]);
        N = Element::mk_n(Element::dof, Element::nnode, Element::n_list[i]);

        detJ = det * Element::w_list[i];

        V += detJ;
        Element::force += (N.row(0)*acc0[0] + N.row(1)*acc0[1]) *detJ;
      }

      Element::force *= Element::mass / V;
    }

  }

// ------------------------------------------------------------------- //
void
  Element::calc_stress() {
    double det;
    Eigen::VectorXd u;
    Eigen::MatrixXd dnj, B;

    std::tie(det, dnj) = Element::mk_dnj(Element::xnT, Element::dn_center);
    B = Element::mk_b(Element::dof, Element::nnode, dnj);
    u = Element::mk_u_hstack(Element::u_p);

    Element::strain = B * u;
    Element::stress = Element::De * Element::strain;
  }

// ------------------------------------------------------------------- //
Eigen::MatrixXd
  Element::mk_m(const Eigen::MatrixXd N) {
    Eigen::MatrixXd M;

    M = N.transpose() * N;
    return M;
  }

Eigen::MatrixXd
  Element::mk_n(const size_t dof, const size_t nnode, const Eigen::VectorXd n) {
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(dof,dof*nnode);

    if (dof == 1) {
      for (size_t i = 0; i < nnode; i++){
        N(0,i) = n(i);
      }
    } else if (dof == 2) {
      for (size_t i = 0; i < nnode; i++){
        size_t i0 = 2*i;
        size_t i1 = 2*i+1;

        N(0,i0) = n(i);
        N(1,i1) = n(i);
      }
    } else if (dof == 3) {
      for (size_t i = 0; i < nnode; i++){
        size_t i0 = 3*i;
        size_t i1 = 3*i+1;
        size_t i2 = 3*i+2;

        N(0,i0) = n(i);
        N(1,i1) = n(i);
        N(2,i2) = n(i);
      }
    }

    return N;
  }

// ------------------------------------------------------------------- //
Eigen::MatrixXd
  Element::mk_nqn(const Eigen::MatrixXd N, const Eigen::MatrixXd q, const Eigen::MatrixXd imp) {
    Eigen::MatrixXd nqn;

    nqn = N.transpose() * q.transpose() * imp * q * N;
    return nqn;
  }

std::tuple<double, Eigen::MatrixXd>
  Element::mk_q(const size_t dof, const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn) {
    Eigen::MatrixXd q;
    Eigen::VectorXd n(2), t(2);
    double det;

    t = xnT * dn;
    n(0) = t(1); n(1) = -t(0);
    det = n.norm();

    if (dof == 1){
      q = Eigen::MatrixXd::Zero(1,1);
      q(0,0) = 1.0;

    } else if (dof == 2){
      q = Eigen::MatrixXd::Zero(2,2);
      q(0,0) = n(0)/det; q(0,1) = n(1)/det;
      q(1,0) = t(0)/det; q(1,1) = t(1)/det;

    } else if (dof == 3){
      q = Eigen::MatrixXd::Zero(3,3);
      q(0,0) = n(0)/det; q(0,1) = n(1)/det;
      q(1,0) = t(0)/det; q(1,1) = t(1)/det;
      q(2,2) = 1.0/det;
    }

    return std::forward_as_tuple(det, q);
  }

// ------------------------------------------------------------------- //
Eigen::MatrixXd
  Element::mk_k(const Eigen::MatrixXd B, const Eigen::MatrixXd D) {
    Eigen::MatrixXd K;

    K = B.transpose() * D * B;
    return K;
  }

Eigen::MatrixXd
  Element::mk_b(const size_t dof, const size_t nnode, const Eigen::MatrixXd dnj) {
    Eigen::MatrixXd B;

    if (dof == 1) {
      B = Eigen::MatrixXd::Zero(2,nnode);
      for (size_t i = 0; i < nnode; i++){
        B(0,i) = dnj(i,0);
        B(1,i) = dnj(i,1);
      }

    } else if (dof == 2) {
      B = Eigen::MatrixXd::Zero(3,2*nnode);
      for (size_t i = 0; i < nnode; i++){
        size_t i0 = 2*i;
        size_t i1 = 2*i+1;

        B(0,i0) = dnj(i,0);
        B(1,i1) = dnj(i,1);
        B(2,i0) = dnj(i,1);
        B(2,i1) = dnj(i,0);
      }

    } else if (dof == 3) {
      B = Eigen::MatrixXd::Zero(5,3*nnode);
      for (size_t i = 0; i < nnode; i++){
        size_t i0 = 3*i;
        size_t i1 = 3*i+1;
        size_t i2 = 3*i+2;

        B(0,i0) = dnj(i,0);
        B(1,i1) = dnj(i,1);
        B(2,i0) = dnj(i,1);
        B(2,i1) = dnj(i,0);

        B(3,i2) = dnj(i,0);
        B(4,i2) = dnj(i,1);
      }
    }

    return B;
  }

Eigen::MatrixXd
  Element::mk_b_T(const size_t dof, const size_t nnode, const Eigen::MatrixXd dnj) {
    Eigen::MatrixXd B;

    if (dof == 1) {
      B = Eigen::MatrixXd::Zero(nnode,2);
      for (size_t i = 0; i < nnode; i++){
        B(i,0) = dnj(i,0);
        B(i,1) = dnj(i,1);
      }

    } else if (dof == 2) {
      B = Eigen::MatrixXd::Zero(2*nnode,3);
      for (size_t i = 0; i < nnode; i++){
        size_t i0 = 2*i;
        size_t i1 = 2*i+1;

        B(i0,0) = dnj(i,0);
        B(i1,1) = dnj(i,1);
        B(i0,2) = dnj(i,1);
        B(i1,2) = dnj(i,0);
      }

    } else if (dof == 3) {
      B = Eigen::MatrixXd::Zero(3*nnode,5);
      for (size_t i = 0; i < nnode; i++){
        size_t i0 = 3*i;
        size_t i1 = 3*i+1;
        size_t i2 = 3*i+2;

        B(i0,0) = dnj(i,0);
        B(i1,1) = dnj(i,1);
        B(i0,2) = dnj(i,1);
        B(i1,2) = dnj(i,0);

        B(i2,3) = dnj(i,0);
        B(i2,4) = dnj(i,1);
      }
    }

    return B;
  }

// ------------------------------------------------------------------- //
Eigen::VectorXd
  Element::Hencky_stress(const Eigen::MatrixXd D, const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    double J;
    Eigen::Matrix2d strain;
    Eigen::VectorXd strain_vector(3), stress;

    std::tie(J, strain) = Euler_log_strain(dnj, u);
    strain_vector << strain(0,0), strain(1,1), strain(0,1)+strain(1,0);
    stress = D * strain_vector / J;
    return stress;
  }

std::tuple<double, Eigen::Matrix2d>
  Element::Euler_log_strain(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    double J;
    Eigen::Matrix2d FF, strain;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es;
    Eigen::Vector2d L, log_L;
    Eigen::Matrix2d P;

    std::tie(J, FF) = Element::mk_FF(dnj, u);
    es.compute(FF);
    L = es.eigenvalues();
    log_L = L.array().log();
    P = es.eigenvectors();

    strain = 0.5 * P * log_L.asDiagonal() * P.transpose();
    return std::forward_as_tuple(J, strain);
  }

std::tuple<double, Eigen::Matrix2d>
  Element::mk_FF(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    double J;
    Eigen::Matrix2d F, FF;

    std::tie(J, F) = Element::mk_F(dnj, u);
    FF = F * F.transpose();
    return std::forward_as_tuple(J, FF);
  }

std::tuple<double, Eigen::Matrix2d>
  Element::mk_F(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    double det;
    Eigen::Matrix2d dnu, F;

    dnu = Element::mk_dnu(dnj, u);
    det = (1.0-dnu(0,0))*(1.0-dnu(1,1)) - dnu(0,1)*dnu(1,0);
    F << 1.0-dnu(1,1),     dnu(0,1),
             dnu(1,0), 1.0-dnu(0,0);
    F /= det;
    return std::forward_as_tuple(1.0/det, F);
  }

Eigen::Matrix2d
  Element::mk_dnu(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    Eigen::Matrix2d dnu;

    dnu = u.transpose() * dnj;
    return dnu;
  }

// ------------------------------------------------------------------- //
std::tuple<double, Eigen::MatrixXd>
  Element::mk_dnj(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn) {
    double det;
    Eigen::Matrix2d jacobi_inv;
    Eigen::MatrixXd dnj;

    std::tie(det, jacobi_inv) = Element::mk_inv_jacobi(xnT, dn);
    dnj = dn * jacobi_inv;
    return std::forward_as_tuple(det, dnj);
  }

std::tuple<double, Eigen::Matrix2d>
  Element::mk_inv_jacobi(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn) {
    double det;
    Eigen::Matrix2d jacobi, jacobi_inv;

    std::tie(det, jacobi) = Element::mk_jacobi(xnT, dn);
    jacobi_inv <<  jacobi(1,1), -jacobi(0,1),
                  -jacobi(1,0),  jacobi(0,0);
    jacobi_inv /= det;
    return std::forward_as_tuple(det, jacobi_inv);
  }

std::tuple<double, Eigen::Matrix2d>
  Element::mk_jacobi(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn) {
    Eigen::Matrix2d jacobi = xnT * dn;

    double det = jacobi(0,0)*jacobi(1,1) - jacobi(0,1)*jacobi(1,0);
    return std::forward_as_tuple(det, jacobi);
  }
