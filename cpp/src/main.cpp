#include "all.h"
#include <Eigen/Core>
#include "element_style.h"
#include "node.h"
#include "element.h"
#include "material.h"
#include "fem.h"
#include "io_data.h"

int main() {

  // ----- Input FEM Mesh ----- //
  Fem fem = input_mesh("input/mesh.in");
  auto outputs = input_outputs("input/output.in");

  // ----- FEM Set up ----- //


}
