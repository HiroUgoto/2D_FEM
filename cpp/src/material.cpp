#include "all.h"
#include "material.h"
#include <Eigen/Core>

Material::Material (int id, std::string style, std::vector<double> param)
  {
    Material::id = id;
    Material::style = style;
    Material::param = param;
  }

void Material::print() {
  std::cout << Material::id << ": ";
  std::cout << Material::style << ", ";
  for (size_t i = 0 ; i < Material::param.size() ; i++) {
    std::cout << Material::param.at(i) << " ";
  }
  std::cout << "\n";
}
