#include "all.h"
#include <Eigen/Dense>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"


// ------------------------------------------------------------------- //
Element::Element (size_t id, std::string style, int material_id, std::vector<size_t> inode)
  {
    this->id = id;
    this->style = style;
    this->material_id = material_id;
    this->inode = inode;

    this->gravity = 9.8;
    this->nnode = inode.size();

    this->set_style();
  }

// ------------------------------------------------------------------- //
void
  Element::print() {
    std::cout << this->id << ": ";
    std::cout << this->style << ", ";
    std::cout << this->material_id << ", ";

    for (size_t i = 0 ; i < this->inode.size() ; i++) {
      std::cout << this->inode.at(i) << " ";
    }
    std::cout << "\n";
  }

void
  Element::set_style() {
    ElementStyle* estyle_p = set_element_style(this->style);

    this->dim = estyle_p->dim;
    this->ng = estyle_p->ng;
    this->xi = estyle_p->xi;
    this->w  = estyle_p->w;

    this->ng_all = estyle_p->ng_all;
    this->n_list  = estyle_p->n_list;
    this->dn_list = estyle_p->dn_list;
    this->w_list = estyle_p->w_list;

    this->dn_center = estyle_p->dn_center;
  }

void
  Element::set_nodes(std::vector<Node*> nodes_p) {
    this->nodes_p = nodes_p;
  }

void
  Element::set_material(Material* material_p) {
    if (material_p == nullptr) {
      this->rho = 0.0;
    } else {
      size_t id = material_p->id;
      std::string style = material_p->style;
      std::vector<double> param = material_p->param;

      this->material.set_init(id,style,param);
      this->rho = this->material.rho;
    }
  }

void
  Element::set_pointer_list(){
    this->u_p.resize(this->nnode);
    this->v_p.resize(this->nnode);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      this->u_p[inode] = &this->nodes_p[inode]->u;
      this->v_p[inode] = &this->nodes_p[inode]->v;
    }
  }

void
  Element::set_xn(){
    this->xnT = Eigen::MatrixXd::Zero(2,this->nnode);

    for (size_t inode = 0 ; inode < this->nnode ; inode++ ) {
      this->xnT(0,inode) = this->nodes_p[inode]->xyz[0] + (*this->u_p[inode])[0];
      this->xnT(1,inode) = this->nodes_p[inode]->xyz[1] + (*this->u_p[inode])[1];
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_local_matrix_init(const size_t dof){
    this->dof = dof;
    this->ndof = dof*this->nnode;

    this->M_diag = Eigen::VectorXd::Zero(this->ndof);

    this->K = Eigen::MatrixXd::Zero(this->ndof,this->ndof);
    this->K_diag = Eigen::VectorXd::Zero(this->ndof);
    this->K_off_diag = Eigen::MatrixXd::Zero(this->ndof,this->ndof);

    this->C = Eigen::MatrixXd::Zero(this->ndof,this->ndof);
    this->C_diag = Eigen::VectorXd::Zero(this->ndof);
    this->C_off_diag = Eigen::MatrixXd::Zero(this->ndof,this->ndof);

    this->force = Eigen::VectorXd::Zero(this->ndof);

    if (this->dim == 2) {
      double V = 0.0;
      for (size_t i = 0 ; i < this->ng_all ; i++){
        auto [det, jacobi] = mk_jacobi(this->xnT, this->dn_list[i]);
        V += det * this->w_list[i];
      }
      this->mass = this->rho * V;

    } else if (this->dim == 1) {
      this->imp = this->material.mk_imp(this->dof);
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_local_matrix() {
    if (this->dim == 2) {
      Eigen::MatrixXd M = Eigen::MatrixXd::Zero(this->ndof,this->ndof);

      this->C = Eigen::MatrixXd::Zero(this->ndof,this->ndof);
      this->K = Eigen::MatrixXd::Zero(this->ndof,this->ndof);

      this->De = this->material.mk_d(this->dof);
      this->Dv = this->material.mk_visco(this->dof);

      for (size_t i = 0 ; i < this->ng_all ; i++){
        double detJ;
        Eigen::MatrixXd N, Me, B, K, C;

        auto [det, dnj] = mk_dnj(this->xnT, this->dn_list[i]);

        N = mk_n(this->dof, this->nnode, this->n_list[i]);
        Me = mk_m(N);

        B = mk_b(this->dof, this->nnode, dnj);
        K = mk_k(B, this->De);
        C = mk_k(B, this->Dv);

        detJ = det * this->w_list[i];

        M += Me * detJ;
        this->C += C * detJ;
        this->K += K * detJ;
      }

      double tr_M = M.trace() / this->dof;
      this->M_diag = M.diagonal() * this->mass/tr_M;

      this->K_diag = this->K.diagonal();
      this->K_off_diag = this->K_diag.asDiagonal();
      this->K_off_diag = (this->K) - (this->K_off_diag);

      this->C_diag = this->C.diagonal();
      this->C_off_diag = this->C_diag.asDiagonal();
      this->C_off_diag = (this->C) - (this->C_off_diag);

    } else if (this->dim == 1) {
      if (this->style.find("input") != std::string::npos) {
        this->C = Eigen::MatrixXd::Zero(this->ndof,this->ndof);

        for (size_t i = 0 ; i < this->ng_all ; i++){
          double detJ;
          Eigen::MatrixXd N, NqN;

          auto [det, q] = mk_q(this->dof, this->xnT,  this->dn_list[i]);

          N = mk_n(this->dof, this->nnode, this->n_list[i]);
          NqN = mk_nqn(N, q, this->imp);

          detJ = det * this->w_list[i];
          this->C += NqN * detJ;
        }

        this->C_diag = this->C.diagonal();
        this->C_off_diag = this->C_diag.asDiagonal();
        this->C_off_diag = (this->C) - (this->C_off_diag);
      }
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_local_vector() {
    if (this->dim == 2) {
      this->force = Eigen::VectorXd::Zero(this->ndof);

      double V = 0.0;
      for (size_t i = 0 ; i < this->ng_all ; i++){
        double detJ;
        Eigen::MatrixXd N;

        auto [det, jacobi] = mk_jacobi(this->xnT, this->dn_list[i]);
        N = mk_n(this->dof, this->nnode, this->n_list[i]);

        detJ = det * this->w_list[i];

        V += detJ;
        this->force += N.row(1)*detJ * this->gravity;
      }

      this->force *= this->mass / V;
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_local_update() {
    if (this->dim == 2) {
      Eigen::MatrixXd M = Eigen::MatrixXd::Zero(this->ndof,this->ndof);
      this->C = Eigen::MatrixXd::Zero(this->ndof,this->ndof);
      this->Dv = this->material.mk_visco(this->dof);
      this->force = Eigen::VectorXd::Zero(this->ndof);

      double V = 0.0;
      for (size_t i = 0 ; i < this->ng_all ; i++){
        double detJ;
        Eigen::MatrixXd N, Me, B, K, C;

        auto [det, dnj] = mk_dnj(this->xnT, this->dn_list[i]);

        N = mk_n(this->dof, this->nnode, this->n_list[i]);
        Me = mk_m(N);

        B = mk_b(this->dof, this->nnode, dnj);
        C = mk_k(B, this->Dv);

        detJ = det * this->w_list[i];
        M += Me * detJ;
        this->C += C * detJ;

        V += detJ;
        this->force += N.row(1)*detJ * this->gravity;
      }

      double tr_M = M.trace() / this->dof;
      this->M_diag = M.diagonal() * this->mass/tr_M;

      this->C_diag = this->C.diagonal();
      this->C_off_diag = this->C_diag.asDiagonal();
      this->C_off_diag = (this->C) - (this->C_off_diag);

      this->force *= this->mass / V;

    } else if (this->dim == 1) {
      if (this->style.find("input") != std::string::npos) {
        this->C = Eigen::MatrixXd::Zero(this->ndof,this->ndof);

        for (size_t i = 0 ; i < this->ng_all ; i++){
          double detJ;
          Eigen::MatrixXd N, NqN;

          auto [det, q] = mk_q(this->dof, this->xnT,  this->dn_list[i]);

          N = mk_n(this->dof, this->nnode, this->n_list[i]);
          NqN = mk_nqn(N, q, this->imp);

          detJ = det * this->w_list[i];
          this->C += NqN * detJ;
        }

        this->C_diag = this->C.diagonal();
        this->C_off_diag = this->C_diag.asDiagonal();
        this->C_off_diag = (this->C) - (this->C_off_diag);
      }
    }
  }


// ------------------------------------------------------------------- //
void
  Element::mk_ku() {
    Eigen::VectorXd u, ku;

    u = this->mk_u_hstack(this->u_p);
    ku = this->K * u;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) += ku(i0+i);
      }
    }
  }

void
  Element::mk_ku_u(const std::vector<Eigen::VectorXd*> u_p) {
    Eigen::VectorXd u, ku;

    u = this->mk_u_hstack(u_p);
    ku = this->K * u;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) += ku(i0+i);
      }
    }

  }

void
  Element::mk_cv() {
    Eigen::VectorXd v, cv;

    v = this->mk_u_hstack(this->v_p);
    cv = this->C_off_diag * v;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) += cv(i0+i);
      }
    }
  }

void
  Element::mk_ku_cv() {
    Eigen::VectorXd u, ku;
    Eigen::VectorXd v, cv;

    u = this->mk_u_hstack(this->u_p);
    v = this->mk_u_hstack(this->v_p);
    ku = this->K * u;
    cv = this->C_off_diag * v;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) += ku(i0+i) + cv(i0+i);
      }
    }
  }

// ------------------------------------------------------------------- //
void
  Element::mk_B_stress() {
    if (this->dim == 1) {
      this->mk_ku();

    } else if (this->dim == 2) {
      Eigen::VectorXd force = Eigen::VectorXd::Zero(this->ndof);
      Eigen::MatrixXd u = this->mk_u_vstack(this->u_p);

      for (size_t i = 0 ; i < this->ng_all ; i++){
        double detJ;
        Eigen::MatrixXd BT;
        Eigen::VectorXd stress;

        auto [det, dnj] = mk_dnj(this->xnT, this->dn_list[i]);
        BT = mk_b_T(this->dof, this->nnode, dnj);
        stress = Hencky_stress(this->De, dnj, u);

        detJ = det * this->w_list[i];
        force += BT * stress * detJ;
      }

      for (size_t inode = 0 ; inode < this->nnode ; inode++){
        size_t i0 = inode*this->dof;
        for (size_t i = 0 ; i < this->dof ; i++) {
          this->nodes_p[inode]->force(i) += force(i0+i);
        }
      }
    }

  }


// ------------------------------------------------------------------- //
Eigen::VectorXd
  Element::mk_u_hstack(const std::vector<Eigen::VectorXd*> u_p) {
    Eigen::VectorXd u(this->ndof);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        u(i0+i) = (*u_p[inode])(i);
      }
    }
    return u;
  }

Eigen::MatrixXd
  Element::mk_u_vstack(const std::vector<Eigen::VectorXd*> u_p) {
    Eigen::MatrixXd u(this->nnode,this->dof);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      for (size_t i = 0 ; i < this->dof ; i++) {
        u(inode,i) = (*u_p[inode])(i);
      }
    }
    return u;
  }


// ------------------------------------------------------------------- //
void
  Element::update_inputwave(const Eigen::VectorXd vel0) {
    Eigen::VectorXd v(this->ndof), cv(this->ndof);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        v(i0+i) = vel0(i);
      }
    }
    cv = this->C * v;

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) -= 2.0*cv(i0+i);
      }
    }
  }

// ------------------------------------------------------------------- //
void
  Element::update_bodyforce(const Eigen::VectorXd acc0) {
    this->mk_bodyforce(acc0);

    for (size_t inode = 0 ; inode < this->nnode ; inode++){
      size_t i0 = inode*this->dof;
      for (size_t i = 0 ; i < this->dof ; i++) {
        this->nodes_p[inode]->force(i) -= this->force(i0+i);
      }
    }
  }

void
  Element::mk_bodyforce(const Eigen::VectorXd acc0) {
    if (this->dim == 2) {
      this->force = Eigen::VectorXd::Zero(this->ndof);

      double V = 0.0;
      for (size_t i = 0 ; i < this->ng_all ; i++){
        double detJ;
        Eigen::MatrixXd N;

        auto [det, jacobi] = mk_jacobi(this->xnT, this->dn_list[i]);
        N = mk_n(this->dof, this->nnode, this->n_list[i]);

        detJ = det * this->w_list[i];

        V += detJ;
        this->force += (N.row(0)*acc0[0] + N.row(1)*acc0[1]) *detJ;
      }

      this->force *= this->mass / V;
    }

  }

// ------------------------------------------------------------------- //
void
  Element::calc_stress() {
    Eigen::VectorXd u;
    Eigen::MatrixXd B;

    auto [det, dnj] = mk_dnj(this->xnT, this->dn_center);
    B = mk_b(this->dof, this->nnode, dnj);
    u = this->mk_u_hstack(this->u_p);

    this->strain = B * u;
    this->stress = this->De * this->strain;
  }


// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
Eigen::MatrixXd
  mk_m(const Eigen::MatrixXd N) {
    Eigen::MatrixXd M;

    M = N.transpose() * N;
    return M;
  }

Eigen::MatrixXd
  mk_n(const size_t dof, const size_t nnode, const Eigen::VectorXd n) {
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
  mk_nqn(const Eigen::MatrixXd N, const Eigen::MatrixXd q, const Eigen::MatrixXd imp) {
    Eigen::MatrixXd nqn;

    nqn = N.transpose() * q.transpose() * imp * q * N;
    return nqn;
  }

std::tuple<double, Eigen::MatrixXd>
  mk_q(const size_t dof, const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn) {
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

    return {det, q};
  }

// ------------------------------------------------------------------- //
Eigen::MatrixXd
  mk_k(const Eigen::MatrixXd B, const Eigen::MatrixXd D) {
    Eigen::MatrixXd K;

    K = B.transpose() * D * B;
    return K;
  }

Eigen::MatrixXd
  mk_b(const size_t dof, const size_t nnode, const Eigen::MatrixXd dnj) {
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
  mk_b_T(const size_t dof, const size_t nnode, const Eigen::MatrixXd dnj) {
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
  Hencky_stress(const Eigen::MatrixXd D, const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    Eigen::VectorXd strain_vector(3), stress;

    auto [J, strain] = Euler_log_strain(dnj, u);
    strain_vector << strain(0,0), strain(1,1), strain(0,1)+strain(1,0);
    stress = D * strain_vector / J;
    return stress;
  }

std::tuple<double, Eigen::Matrix2d>
  Euler_log_strain(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    Eigen::Matrix2d strain;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es;
    Eigen::Vector2d L, log_L;
    Eigen::Matrix2d P;

    auto [J, FF] = mk_FF(dnj, u);
    es.compute(FF);
    L = es.eigenvalues();
    log_L = L.array().log();
    P = es.eigenvectors();

    strain = 0.5 * P * log_L.asDiagonal() * P.transpose();
    return {J, strain};
  }

std::tuple<double, Eigen::Matrix2d>
  mk_FF(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    Eigen::Matrix2d FF;

    auto [J, F] = mk_F(dnj, u);
    FF = F * F.transpose();
    return std::forward_as_tuple(J, FF);
  }

std::tuple<double, Eigen::Matrix2d>
  mk_F(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    double det;
    Eigen::Matrix2d dnu, F;

    dnu = mk_dnu(dnj, u);
    det = (1.0-dnu(0,0))*(1.0-dnu(1,1)) - dnu(0,1)*dnu(1,0);
    F << 1.0-dnu(1,1),     dnu(0,1),
             dnu(1,0), 1.0-dnu(0,0);
    F /= det;
    return {1.0/det, F};
  }

Eigen::Matrix2d
  mk_dnu(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u) {
    Eigen::Matrix2d dnu;

    dnu = u.transpose() * dnj;
    return dnu;
  }

// ------------------------------------------------------------------- //
std::tuple<double, Eigen::MatrixXd>
  mk_dnj(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn) {
    Eigen::MatrixXd dnj;

    auto [det, jacobi_inv] = mk_inv_jacobi(xnT, dn);
    dnj = dn * jacobi_inv;
    return {det, dnj};
  }

std::tuple<double, Eigen::Matrix2d>
  mk_inv_jacobi(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn) {
    Eigen::Matrix2d jacobi_inv;

    auto [det, jacobi] = mk_jacobi(xnT, dn);
    jacobi_inv <<  jacobi(1,1), -jacobi(0,1),
                  -jacobi(1,0),  jacobi(0,0);
    jacobi_inv /= det;
    return {det, jacobi_inv};
  }


std::tuple<double, Eigen::Matrix2d>
  mk_jacobi(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn) {
    Eigen::Matrix2d jacobi = xnT * dn;

    double det = jacobi(0,0)*jacobi(1,1) - jacobi(0,1)*jacobi(1,0);
    return {det, jacobi};
    // return std::forward_as_tuple(det, jacobi);
  }
