#include "all.h"
#include <Eigen/Core>
#include "element_style.h"
#include "element.h"


Element::Element (size_t id, std::string style, size_t material_id, std::vector<size_t> inode)
  {
    Element::id = id;
    Element::style = style;
    Element::material_id = material_id;
    Element::inode = inode;

    Element::gravity = 9.8;
    Element::nnode = inode.size();

    Element::set_style();
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

void Element::set_style() {
  Element::estyle_p = set_element_style(Element::style);
  Element::dim = Element::estyle_p->dim;
  
  Element::ng = Element::estyle_p->ng;
  Element::xi = Element::estyle_p->xi;
  Element::w  = Element::estyle_p->w;
}
