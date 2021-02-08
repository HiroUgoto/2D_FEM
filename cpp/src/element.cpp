#include "all.h"
#include "element.h"
#include <Eigen/Core>


Element::Element (int id, std::string style, int material_id, std::vector<int> inode)
  {
    Element::id = id;
    Element::style = style;
    Element::material_id = material_id;
    Element::inode = inode;
  }

void Element::print() {
  std::cout << Element::id << ": ";
  std::cout << Element::style << ", ";
  std::cout << Element::material_id << ", ";

  for (size_t i = 0 ; i < Element::inode.size() ; i++) {
    std::cout << Element::inode.at(i) << " ";
  }
  std::cout << "\n";
}
