#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "fem.h"

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
Fem::Fem (size_t dof, std::vector<Node> nodes,
                std::vector<Element> elements,
                std::vector<Material> materials)
  {
    Fem::nnode = nodes.size();
    Fem::nelem = elements.size();
    Fem::dof = dof;

    Fem::nodes = nodes;
    Fem::elements = elements;
    Fem::materials = materials;

  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::set_init() {
  Fem::_set_mesh();
  Fem::_set_initial_condition();
  Fem::_set_initial_matrix();
}

// ----------------------------- //
void Fem::_set_mesh() {
  for (size_t ielem = 0 ; ielem < Fem::nelem ; ielem++ ){
    Element& element = Fem::elements[ielem];

    std::vector<Node*> nodes_p(element.nnode);
    for (size_t i = 0 ; i < element.nnode ; i++ ){
      size_t inode = element.inode[i];

      if (Fem::nodes[inode].id == inode) {
        nodes_p[i] = &Fem::nodes[inode];
      } else {
        for (size_t i0 = 0 ; i0 < Fem::nnode ; i0++ ) {
          if (Fem::nodes[i0].id == inode) {
            nodes_p[i] = &Fem::nodes[i0];
            break;
          }
        }
      }

    }

    element.set_nodes(nodes_p);

    Material* material_p = nullptr;
    if (Fem::materials[element.material_id].id == element.material_id) {
      material_p = &Fem::materials[element.material_id];

    } else {
      for (size_t i0 = 0 ; i0 < Fem::materials.size() ; i0++) {
        if (Fem::materials[i0].id == element.material_id) {
          material_p = &Fem::materials[i0];
          break;
        }
      }
    }

    element.set_material(material_p);


    if (element.style.find("input") != std::string::npos) {
      Fem::input_elements.push_back(ielem);
    }
    if (element.style.find("connect") != std::string::npos) {
      Fem::connected_elements.push_back(ielem);
    }

  }

}

// ----------------------------- //
void Fem::_set_initial_condition() {
  for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
    Node& node = Fem::nodes[inode];

    node.set_initial_condition();

    size_t total_dof = std::accumulate(node.freedom.begin(),node.freedom.end(),0);
    if (total_dof == dof) {
      Fem::free_nodes.push_back(inode);
    } else {
      Fem::fixed_nodes.push_back(inode);
    }
  }

  for (size_t ielem = 0 ; ielem < Fem::nelem ; ielem++ ){
    Element& element = Fem::elements[ielem];
    element.set_pointer_list();
  }
}

// ----------------------------- //
void Fem::_set_initial_matrix(){
  for (size_t ielem = 0 ; ielem < Fem::nelem ; ielem++ ){
    Element& element = Fem::elements[ielem];

    element.set_xn();
    element.mk_local_matrix_init(Fem::dof);
    element.mk_local_matrix();
    element.mk_local_vector();

    size_t id = 0;
    for (size_t inode = 0 ; inode < element.nnode ; inode++) {
      for (size_t i = 0 ; i < Fem::dof ; i++) {
        element.nodes_p[inode]->mass[i] += element.M_diag[id];
        element.nodes_p[inode]->c[i] += element.C_diag[id];
        element.nodes_p[inode]->k[i] += element.K_diag[id];
        element.nodes_p[inode]->static_force[i] += element.force[id];
        id++;
      }
    }
  }
}

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::set_output(std::tuple<std::vector<size_t>, std::vector<size_t>> outputs) {
  std::vector<size_t> output_node_list, output_element_list;
  std::tie(output_node_list, output_element_list) = outputs;

  Fem::output_nnode = output_node_list.size();
  Fem::output_nodes = output_node_list;

  Fem::output_nelem = output_element_list.size();
  Fem::output_elements = output_element_list;
}
