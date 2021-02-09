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
