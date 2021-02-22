#include "all.h"
#include <Eigen/Core>
#include "element_style.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;


// ------- Set element style ---------------------------- //
ElementStyle* set_element_style(const std::string style) {
  ElementStyle* es_p = nullptr;

  if (style == "2d4solid") {
    es_p = new Solid_2d_4Node();
  } else if (style == "2d9solid") {
    es_p = new Solid_2d_9Node();
  } else if (style == "1d2line") {
    es_p = new Line_1d_2Node();
  } else if (style == "1d3line") {
    es_p = new Line_1d_3Node();
  } else if (style == "1d2input") {
    es_p = new Input_1d_2Node();
  } else if (style == "1d3input") {
    es_p = new Input_1d_3Node();
  } else if (style == "connect") {
    es_p = new Connect();
  } else {
    es_p = new ElementStyle();
  }
  return es_p;
}

// ------- Set gauss points ---------------------------- //
void set_gauss_points (const size_t n, EV& xi, EV& w) {
  if (n == 3){
    xi = EV::Zero(3);
    w  = EV::Zero(3);

    xi << 0.0, sqrt(15.0)/5.0, -sqrt(15.0)/5.0;
    w  << 8.0/9.0, 5.0/9.0, 5.0/9.0;

  } else if (n == 5) {
    xi = EV::Zero(5);
    w  = EV::Zero(5);

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

EV ElementStyle::shape_function_n (double xi, double zeta) {
  EV n;
  return n;
}

EM ElementStyle::shape_function_dn (double xi, double zeta) {
  EM nd;
  return nd;
}

// ------- Actual element style class ------------------- //
// ----------------------------------------------------- //
Solid_2d_4Node::Solid_2d_4Node () {
  this->dim = 2;
  this->ng = 3;
  set_gauss_points(this->ng, this->xi, this->w);

  this->ng_all = this->ng * this->ng;
  this->n_list.resize(this->ng_all);
  this->dn_list.resize(this->ng_all);
  this->w_list.resize(this->ng_all);

  size_t id = 0;
  for (size_t i=0 ; i<ng ; i++) {
    for (size_t j=0 ; j<ng ; j++) {
      n_list[id] = this->shape_function_n(this->xi[i],this->xi[j]);
      dn_list[id] = this->shape_function_dn(this->xi[i],this->xi[j]);
      w_list[id] = w[i]*w[j];
      id++;
    }
  }

  dn_center = this->shape_function_dn(0.0,0.0);
}

EV Solid_2d_4Node::shape_function_n (double xi, double zeta) {
    EV n = EV::Zero(4);
    n(0) = (1.0 - xi)*(1.0 - zeta) / 4.0;
    n(1) = (1.0 + xi)*(1.0 - zeta) / 4.0;
    n(2) = (1.0 + xi)*(1.0 + zeta) / 4.0;
    n(3) = (1.0 - xi)*(1.0 + zeta) / 4.0;
    return n;
  }

EM Solid_2d_4Node::shape_function_dn (double xi, double zeta) {
    EM dn = EM::Zero(4,2);
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
Solid_2d_8Node::Solid_2d_8Node () {
  this->dim = 2;
  this->ng = 5;
  set_gauss_points(this->ng, this->xi, this->w);

  this->ng_all = this->ng * this->ng;
  this->n_list.resize(this->ng_all);
  this->dn_list.resize(this->ng_all);
  this->w_list.resize(this->ng_all);

  size_t id = 0;
  for (size_t i=0 ; i<ng ; i++) {
    for (size_t j=0 ; j<ng ; j++) {
      n_list[id] = this->shape_function_n(this->xi[i],this->xi[j]);
      dn_list[id] = this->shape_function_dn(this->xi[i],this->xi[j]);
      w_list[id] = w[i]*w[j];
      id++;
    }
  }

  dn_center = this->shape_function_dn(0.0,0.0);
}

EV Solid_2d_8Node::shape_function_n (double xi, double zeta) {
    EV n = EV::Zero(8);
    n(0) = (1.0 - xi)*(1.0 - zeta)*(-1.0-xi-zeta) / 4.0;
    n(1) = (1.0 + xi)*(1.0 - zeta)*(-1.0+xi-zeta) / 4.0;
    n(2) = (1.0 + xi)*(1.0 + zeta)*(-1.0+xi+zeta) / 4.0;
    n(3) = (1.0 - xi)*(1.0 + zeta)*(-1.0-xi+zeta) / 4.0;

    n(4) = (1.0 - xi*xi)*(1.0 - zeta) / 2.0;
    n(5) = (1.0 + xi)*(1.0 - zeta*zeta) / 2.0;
    n(6) = (1.0 - xi*xi)*(1.0 + zeta) / 2.0;
    n(7) = (1.0 - xi)*(1.0 - zeta*zeta) / 2.0;
    return n;
  }

EM Solid_2d_8Node::shape_function_dn (double xi, double zeta) {
    EM dn = EM::Zero(8,2);
    dn(0,0) = (1.0 - zeta)*(2.0*xi+zeta) / 4.0;
    dn(0,1) = (1.0 -   xi)*(xi+2.0*zeta) / 4.0;

    dn(1,0) = (1.0 - zeta)*(2.0*xi-zeta) / 4.0;
    dn(1,1) = -(1.0 +  xi)*(xi-2.0*zeta) / 4.0;

    dn(2,0) = (1.0 + zeta)*(2.0*xi+zeta) / 4.0;
    dn(2,1) = (1.0 +   xi)*(xi+2.0*zeta) / 4.0;

    dn(3,0) = (1.0 + zeta)*(2.0*xi-zeta) / 4.0;
    dn(3,1) = -(1.0 -  xi)*(xi-2.0*zeta) / 4.0;

    dn(4,0) = -xi*(1.0 - zeta);
    dn(4,1) = (xi*xi-1.0) / 2.0;

    dn(5,0) = (1.0 - zeta*zeta) / 2.0;
    dn(5,1) = -(1.0 + xi)*zeta;

    dn(6,0) = -xi*(1.0 + zeta);
    dn(6,1) = (1.0 - xi*xi) / 2.0;

    dn(7,0) = -(1.0 - zeta*zeta) / 2.0;
    dn(7,1) = -(1.0 - xi)*zeta;
    return dn;
  }

// ----------------------------------------------------- //
Solid_2d_9Node::Solid_2d_9Node () {
  this->dim = 2;
  this->ng = 5;
  set_gauss_points(this->ng, this->xi, this->w);

  this->ng_all = this->ng * this->ng;
  this->n_list.resize(this->ng_all);
  this->dn_list.resize(this->ng_all);
  this->w_list.resize(this->ng_all);

  size_t id = 0;
  for (size_t i=0 ; i<ng ; i++) {
    for (size_t j=0 ; j<ng ; j++) {
      n_list[id] = this->shape_function_n(this->xi[i],this->xi[j]);
      dn_list[id] = this->shape_function_dn(this->xi[i],this->xi[j]);
      w_list[id] = w[i]*w[j];
      id++;
    }
  }

  dn_center = this->shape_function_dn(0.0,0.0);
}

EV Solid_2d_9Node::shape_function_n (double xi, double zeta) {
    EV n = EV::Zero(9);
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

EM Solid_2d_9Node::shape_function_dn (double xi, double zeta) {
    EM dn = EM::Zero(9,2);
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

// ----------------------------------------------------- //
Line_1d_2Node::Line_1d_2Node () {
  this->dim = 1;
  this->ng = 3;
  set_gauss_points(this->ng, this->xi, this->w);

  this->ng_all = this->ng;
  this->n_list.resize(this->ng_all);
  this->dn_list.resize(this->ng_all);
  this->w_list.resize(this->ng_all);

  size_t id = 0;
  for (size_t i=0 ; i<ng ; i++) {
    n_list[id] = this->shape_function_n(this->xi[i]);
    dn_list[id] = this->shape_function_dn(this->xi[i]);
    w_list[id] = w[i];
    id++;
  }

  dn_center = this->shape_function_dn(0.0);
}

EV Line_1d_2Node::shape_function_n (double xi, double zeta) {
    EV n = EV::Zero(2);
    n(0) = (1.0 - xi) / 2.0;
    n(1) = (1.0 + xi) / 2.0;
    return n;
  }

EM Line_1d_2Node::shape_function_dn (double xi, double zeta) {
    EM dn = EM::Zero(2,1);
    dn(0,0) = -0.5;
    dn(1,0) =  0.5;
    return dn;
  }

// ----------------------------------------------------- //
Line_1d_3Node::Line_1d_3Node () {
  this->dim = 1;
  this->ng = 5;
  set_gauss_points(this->ng, this->xi, this->w);

  this->ng_all = this->ng;
  this->n_list.resize(this->ng_all);
  this->dn_list.resize(this->ng_all);
  this->w_list.resize(this->ng_all);

  size_t id = 0;
  for (size_t i=0 ; i<ng ; i++) {
    n_list[id] = this->shape_function_n(this->xi[i]);
    dn_list[id] = this->shape_function_dn(this->xi[i]);
    w_list[id] = w[i];
    id++;
  }

  dn_center = this->shape_function_dn(0.0);
}

EV Line_1d_3Node::shape_function_n (double xi, double zeta) {
    EV n = EV::Zero(3);
    n(0) = -xi*(1.0 - xi) / 2.0;
    n(1) =  xi*(1.0 + xi) / 2.0;
    n(2) = (1.0 - xi)*(1.0 + xi);
    return n;
  }

EM Line_1d_3Node::shape_function_dn (double xi, double zeta) {
    EM dn = EM::Zero(3,1);
    dn(0,0) = xi - 0.5;
    dn(1,0) = xi + 0.5;
    dn(2,0) = -2.0*xi;
    return dn;
  }

// ----------------------------------------------------- //
Connect::Connect () {
  this->dim = 0;
  this->ng = 1;
  this->xi = EV::Zero(1);
  this->w  = EV::Zero(1);
}
