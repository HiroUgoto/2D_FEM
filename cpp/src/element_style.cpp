#include "all.h"
#include <Eigen/Core>
#include "element_style.h"

// ------- Set element style ---------------------------- //
ElementStyle* set_element_style(const std::string style) {

  if (style == "2d4solid") {
    ElementStyle* es_p = new Solid_2d_4Node();
    return es_p;

  } else if (style == "2d9solid") {
    ElementStyle* es_p = new Solid_2d_9Node();
    return es_p;

  } else {
    ElementStyle* es_p = new ElementStyle();
    return es_p;
  }

}

// ------- Set gauss points ---------------------------- //
void set_gauss_points (const size_t n, Eigen::VectorXd& xi, Eigen::VectorXd& w) {
  if (n == 3){
    xi = Eigen::VectorXd::Zero(3);
    w  = Eigen::VectorXd::Zero(3);

    xi << 0.0, sqrt(15.0)/5.0, -sqrt(15.0)/5.0;
    w  << 8.0/9.0, 5.0/9.0, 5.0/9.0;

  } else if (n == 5) {
    xi = Eigen::VectorXd::Zero(5);
    w  = Eigen::VectorXd::Zero(5);

    xi << 0.0,
          sqrt(245.0 - 14.0*sqrt(70.0))/21.0,
         -sqrt(245.0 - 14.0*sqrt(70.0))/21.0,
          sqrt(245.0 + 14.0*sqrt(70.0))/21.0,
         -sqrt(245.0 + 14.0*sqrt(70.0))/21.0;
    w << 128.0/225.0,
        (322.0 + 13.0*sqrt(70.0))/900.0,
        (322.0 + 13.0*sqrt(70.0))/900.0,
        (322.0 - 13.0*sqrt(70.0))/900.0,
        (322.0 - 13.0*sqrt(70.0))/900.0;
  }

}




// ------- Virtual class (ElementStyle) ----------------- //
ElementStyle::ElementStyle () {}

Eigen::VectorXd ElementStyle::shape_function_n (double xi, double zeta) {
  Eigen::VectorXd n;
  return n;
}

Eigen::MatrixXd ElementStyle::shape_function_dn (double xi, double zeta) {
  Eigen::MatrixXd nd;
  return nd;
}

// ------- Actual element style class ------------------- //
// ----------------------------------------------------- //
Solid_2d_4Node::Solid_2d_4Node () {
  Solid_2d_4Node::dim = 2;

  size_t ng = 3;
  Solid_2d_4Node::ng = ng;
  set_gauss_points(Solid_2d_4Node::ng, Solid_2d_4Node::xi, Solid_2d_4Node::w);

  Solid_2d_4Node::ng_all = ng*ng;
  Solid_2d_4Node::n_list.resize(ng*ng);
  Solid_2d_4Node::dn_list.resize(ng*ng);
  Solid_2d_4Node::w_list.resize(ng*ng);

  size_t id = 0;
  for (size_t i=0 ; i<ng ; i++) {
    for (size_t j=0 ; j<ng ; j++) {
      Solid_2d_4Node::n_list[id] = Solid_2d_4Node::shape_function_n(Solid_2d_4Node::xi[i],Solid_2d_4Node::xi[j]);
      Solid_2d_4Node::dn_list[id] = Solid_2d_4Node::shape_function_dn(Solid_2d_4Node::xi[i],Solid_2d_4Node::xi[j]);
      Solid_2d_4Node::w_list[id] = Solid_2d_4Node::w[i]*Solid_2d_4Node::w[j];
      id++;
    }
  }
}

Eigen::VectorXd Solid_2d_4Node::shape_function_n (double xi, double zeta) {
  Eigen::VectorXd n = Eigen::VectorXd::Zero(4);
  n(0) = (1.0 - xi)*(1.0 - zeta) / 4.0;
  n(1) = (1.0 + xi)*(1.0 - zeta) / 4.0;
  n(2) = (1.0 + xi)*(1.0 + zeta) / 4.0;
  n(3) = (1.0 - xi)*(1.0 + zeta) / 4.0;
  return n;
}

Eigen::MatrixXd Solid_2d_4Node::shape_function_dn (double xi, double zeta) {
  Eigen::MatrixXd dn = Eigen::MatrixXd::Zero(4,2);
  dn(0,0) = -(1.0 - zeta) / 4.0;
  dn(0,1) = -(1.0 -   xi) / 4.0;

  dn(1,0) =  (1.0 - zeta) / 4.0;
  dn(1,1) = -(1.0 +   xi) / 4.0;

  dn(2,0) =  (1.0 + zeta) / 4.0;
  dn(2,1) =  (1.0 +   xi) / 4.0;

  dn(3,0) = -(1.0 + zeta) / 4.0;
  dn(3,1) =  (1.0 -   xi) / 4.0;
  return dn;
}


// ----------------------------------------------------- //
Solid_2d_9Node::Solid_2d_9Node () {
  Solid_2d_9Node::dim = 2;

  size_t ng = 5;
  Solid_2d_9Node::ng = ng;
  set_gauss_points(Solid_2d_9Node::ng, Solid_2d_9Node::xi, Solid_2d_9Node::w);

  Solid_2d_9Node::ng_all = ng*ng;
  Solid_2d_9Node::n_list.resize(ng*ng);
  Solid_2d_9Node::dn_list.resize(ng*ng);
  Solid_2d_9Node::w_list.resize(ng*ng);

  size_t id = 0;
  for (size_t i=0 ; i<ng ; i++) {
    for (size_t j=0 ; j<ng ; j++) {
      Solid_2d_9Node::n_list[id] = Solid_2d_9Node::shape_function_n(Solid_2d_9Node::xi[i],Solid_2d_9Node::xi[j]);
      Solid_2d_9Node::dn_list[id] = Solid_2d_9Node::shape_function_dn(Solid_2d_9Node::xi[i],Solid_2d_9Node::xi[j]);
      Solid_2d_9Node::w_list[id] = Solid_2d_9Node::w[i]*Solid_2d_9Node::w[j];
      id++;
    }
  }
}

Eigen::VectorXd Solid_2d_9Node::shape_function_n (double xi, double zeta) {
  Eigen::VectorXd n = Eigen::VectorXd::Zero(9);
  n(0) =  (1.0 - xi)*(1.0 - zeta)*xi*zeta / 4.0;
  n(1) = -(1.0 + xi)*(1.0 - zeta)*xi*zeta / 4.0;
  n(2) =  (1.0 + xi)*(1.0 + zeta)*xi*zeta / 4.0;
  n(3) = -(1.0 - xi)*(1.0 + zeta)*xi*zeta / 4.0;

  n(4) = -(1.0 - xi*xi)*zeta*(1.0 - zeta) / 2.0;
  n(5) =  (1.0 + xi)*xi*(1.0 - zeta*zeta) / 2.0;
  n(6) =  (1.0 - xi*xi)*zeta*(1.0 + zeta) / 2.0;
  n(7) = -(1.0 - xi)*xi*(1.0 - zeta*zeta) / 2.0;

  n(8) = (1.0-xi*xi)*(1.0-zeta*zeta);
  return n;
}

Eigen::MatrixXd Solid_2d_9Node::shape_function_dn (double xi, double zeta) {
  Eigen::MatrixXd dn = Eigen::MatrixXd::Zero(9,2);
  dn(0,0) =  (2.0*xi-1.0)*(zeta-1.0)*zeta / 4.0;
  dn(0,1) =  (xi-1.0)*xi*(2.0*zeta-1.0) / 4.0;

  dn(1,0) =  (2.0*xi+1.0)*(zeta-1.0)*zeta / 4.0;
  dn(1,1) =  (xi+1.0)*xi*(2.0*zeta-1.0) / 4.0;

  dn(2,0) =  (2.0*xi+1.0)*(zeta+1.0)*zeta / 4.0;
  dn(2,1) =  (xi+1.0)*xi*(2.0*zeta+1.0) / 4.0;

  dn(3,0) =  (2.0*xi-1.0)*(zeta+1.0)*zeta / 4.0;
  dn(3,1) =  (xi-1.0)*xi*(2.0*zeta+1.0) / 4.0;

  dn(4,0) =  xi*(1.0-zeta)*zeta;
  dn(4,1) =  (1.0-xi*xi)*(2.0*zeta-1.0) / 2.0;

  dn(5,0) =  (2.0*xi+1.0)*(1.0-zeta*zeta) / 2.0;
  dn(5,1) =  -xi*(1.0+xi)*zeta;

  dn(6,0) =  -xi*(1.0+zeta)*zeta;
  dn(6,1) =  (1.0-xi*xi)*(2.0*zeta+1.0) / 2.0;

  dn(7,0) =  (2.0*xi-1.0)*(1.0-zeta*zeta) / 2.0;
  dn(7,1) =  xi*(1.0-xi)*zeta;

  dn(8,0) =  -2.0*xi*(1.0-zeta*zeta);
  dn(8,1) =  -2.0*(1.0-xi*xi)*zeta;
  return dn;
}
