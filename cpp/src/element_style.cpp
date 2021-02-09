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

// ------- Actual element style class ------------------- //
// ----------------------------------------------------- //
Solid_2d_4Node::Solid_2d_4Node () {
  Solid_2d_4Node::dim = 2;
  Solid_2d_4Node::ng = 3;
  set_gauss_points(Solid_2d_4Node::ng, Solid_2d_4Node::xi, Solid_2d_4Node::w);
}

Eigen::VectorXd Solid_2d_4Node::shape_function_n (double xi, double zeta) {
  Eigen::VectorXd n(4);
  n(0) = (1.0 - xi)*(1.0 - zeta) / 4.0;
  n(1) = (1.0 + xi)*(1.0 - zeta) / 4.0;
  n(2) = (1.0 + xi)*(1.0 + zeta) / 4.0;
  n(3) = (1.0 - xi)*(1.0 + zeta) / 4.0;
  return n;
}


// ----------------------------------------------------- //
Solid_2d_9Node::Solid_2d_9Node () {
  Solid_2d_9Node::dim = 2;
  Solid_2d_9Node::ng = 5;
  set_gauss_points(Solid_2d_9Node::ng, Solid_2d_9Node::xi, Solid_2d_9Node::w);
}

Eigen::VectorXd Solid_2d_9Node::shape_function_n (double xi, double zeta) {
  Eigen::VectorXd n(9);
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
