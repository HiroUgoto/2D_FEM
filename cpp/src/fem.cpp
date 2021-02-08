#include "all.h"
#include "node.h"
#include "element.h"
#include "material.h"
#include "fem.h"
#include <Eigen/Core>

Fem::Fem (int dof, std::vector<Node> nodes,
              std::vector<Element> elements,
              std::vector<Material> materials)
  {
    Fem::dof = dof;
    Fem::nodes = nodes;
    Fem::elements = elements;
    Fem::materials = materials;
  }
