#include "all.h"
#include <Eigen/Core>
#include "element_style.h"
#include "node.h"
#include "element.h"
#include "material.h"
#include "fem.h"

Fem::Fem (size_t dof, std::vector<Node> nodes,
                std::vector<Element> elements,
                std::vector<Material> materials)
  {
    Fem::dof = dof;
    Fem::nodes = nodes;
    Fem::elements = elements;
    Fem::materials = materials;
  }
