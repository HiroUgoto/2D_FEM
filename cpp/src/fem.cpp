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
void
  Fem::set_init() {
    Fem::_set_mesh();
    Fem::_set_initial_condition();
    Fem::_set_initial_matrix();
  }

// ----------------------------- //
void
  Fem::_set_mesh() {
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
void
  Fem::_set_initial_condition() {
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
void
  Fem::_set_initial_matrix(){
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
void
  Fem::set_output(std::tuple<std::vector<size_t>, std::vector<size_t>> outputs) {
    std::vector<size_t> output_node_list, output_element_list;
    std::tie(output_node_list, output_element_list) = outputs;

    Fem::output_nnode = output_node_list.size();
    Fem::output_nodes = output_node_list;

    Fem::output_nelem = output_element_list.size();
    Fem::output_elements = output_element_list;
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void
  Fem::self_gravity() {
    // Initial condition //
    double H = 0.0;
    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      if (Fem::nodes[inode].xyz[1] > H) {
        H = Fem::nodes[inode].xyz[1];
      }
    }

    double g = 9.8;
    double vp = 1500.0;
    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      Node& node = Fem::nodes[inode];

      node.u[0] = 0.0;
      node.u[1] = g/(2*vp*vp) * (H*H - node.xyz[1]*node.xyz[1]);
      node.um = node.u;
    }

    Fem::_self_gravity_cg(false);
    Fem::_self_gravity_cg(true);

    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      Node& node = Fem::nodes[inode];
      node.u0 = node.u;
    }

  }

// ------------------------------------------------------------------- //
void
  Fem::_self_gravity_cg(const bool full) {
    size_t id;
    if (full) {
      id = 0;
    } else {
      id = 1;
    }

    // CG Method //
    // --- set initial variables --- //
    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      Node& node = Fem::nodes[inode];
      node.force = Eigen::VectorXd::Zero(node.dof);
    }
    for (size_t ielem = 0 ; ielem < Fem::nelem ; ielem++) {
      Element& element = Fem::elements[ielem];
      element.mk_ku();
    }

    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      Node& node = Fem::nodes[inode];
      for (size_t i = id ; i < node.dof ; i++) {
        if (node.freedom[i] == 0) {
          node._ur[i] = 0.0;
        } else {
          node._ur[i] = node.static_force[i] - node.force[i];
        }
      }
    }
    for (size_t i = 0 ; i < Fem::connected_elements.size() ; i++) {
      size_t ielem = Fem::connected_elements[i];
      Element& element = Fem::elements[ielem];
      Eigen::VectorXd u = Eigen::VectorXd::Zero(element.dof);

      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        u += element.nodes_p[inode]->_ur;
      }
      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        element.nodes_p[inode]->_ur = u/element.nnode;
      }
    }

    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      Node& node = Fem::nodes[inode];
      node._up = node._ur;
    }

    for (size_t ielem = 0 ; ielem < Fem::nelem ; ielem++) {
      Element& element = Fem::elements[ielem];

      element._up_p.resize(element.nnode);
      for (size_t inode = 0 ; inode < element.nnode ; inode++){
        element._up_p[inode] = &element.nodes_p[inode]->_up;
      }
    }

    // --- CG iterations --- //
    for (size_t it=0 ; it < 10*Fem::nnode ; it++) {
      double rr, rr1, py, alpha, beta;

      // y = Ap
      for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
        Node& node = Fem::nodes[inode];
        node.force = Eigen::VectorXd::Zero(node.dof);
      }
      for (size_t ielem = 0 ; ielem < Fem::nelem ; ielem++) {
        Element& element = Fem::elements[ielem];
        element.mk_ku_u(element._up_p);
      }
      for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
        Node& node = Fem::nodes[inode];
        node._uy = node.force;
      }

      // correction boundary condition
      for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
        Node& node = Fem::nodes[inode];
        for (size_t i = id ; i < node.dof ; i++) {
          if (node.freedom[i] == 0) {
            node._uy[i] = 0.0;
          }
        }
      }
      for (size_t i = 0 ; i < Fem::connected_elements.size() ; i++) {
        size_t ielem = Fem::connected_elements[i];
        Element& element = Fem::elements[ielem];
        Eigen::VectorXd u = Eigen::VectorXd::Zero(element.dof);

        for (size_t inode = 0 ; inode < element.nnode ; inode++) {
          u += element.nodes_p[inode]->_uy;
        }
        for (size_t inode = 0 ; inode < element.nnode ; inode++) {
          element.nodes_p[inode]->_uy = u/element.nnode;
        }
      }

      // alpha = rr/py
      rr = 0.0; py = 0.0;
      for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
        Node& node = Fem::nodes[inode];
        rr += node._ur.dot(node._ur);
        py += node._up.dot(node._uy);
      }
      alpha = rr/py;

      // x = x + alpha*p
      rr1 = 0.0;
      for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
        Node& node = Fem::nodes[inode];
        for (size_t i = id ; i < node.dof ; i++) {
          if (node.freedom[i] != 0) {
            node.u[i] += alpha * node._up[i];
            node._ur[i] -= alpha * node._uy[i];
          }
        }
        rr1 += node._ur.dot(node._ur);
      }

      if (rr1 < 1.e-10) {
        break;
      }

      // p = r + beta*p
      beta = rr1/rr;
      for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
        Node& node = Fem::nodes[inode];
        for (size_t i = id ; i < node.dof ; i++) {
          if (node.freedom[i] != 0) {
            node._up[i] = node._ur[i] + beta * node._up[i];
          }
        }
      }

      if (it%100 == 0){
        std::cout << " (self gravity process .. ) " << it << " " ;
        std::cout << Fem::nodes[0].u[1] << " ";
        std::cout << rr1 << "\n";
      }
    }

  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void
  Fem::update_init(const double dt) {
    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      Node& node = Fem::nodes[inode];

      for (size_t i = 0 ; i < node.dof ; i++) {
        node.inv_mc[i] = 1.0 / (node.mass[i] + 0.5*dt*node.c[i]);
        node.mass_inv_mc[i] = node.mass[i] * node.inv_mc[i];
        node.c_inv_mc[i] = node.c[i] * node.inv_mc[i] * 0.5*dt;
        node.dtdt_inv_mc[i] = dt*dt*node.inv_mc[i];
      }
    }

    Fem::dt = dt;
    Fem::inv_dt2 = 1.0/(2.0*dt);
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void
  Fem::update_time(const Eigen::VectorXd acc0) {
    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      Node& node = Fem::nodes[inode];
      node.dynamic_force = node.static_force;
    }

    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      Node& node = Fem::nodes[inode];
      node.force = -node.dynamic_force;
    }

    for (size_t ielem = 0 ; ielem < Fem::nelem ; ielem++) {
      Element& element = Fem::elements[ielem];
      element.update_bodyforce(acc0);
    }

    for (size_t ielem = 0 ; ielem < Fem::nelem ; ielem++) {
      Element& element = Fem::elements[ielem];
      element.mk_ku_cv();
    }

    for (size_t inode = 0 ; inode < Fem::nnode ; inode++) {
      Node& node = Fem::nodes[inode];
      Eigen::VectorXd u = node.u;

      for (size_t i = 0 ; i < node.dof ; i++) {
        if (node.freedom[i] == 0) {
          node.u[i] = 0.0;
        } else {
          node.u[i] = node.mass_inv_mc[i]*(2.0*u[i]-node.um[i]) + node.c_inv_mc[i]*node.um[i] - node.dtdt_inv_mc[i]*node.force[i];
        }
      }
      node.v = (node.u - node.um) * Fem::inv_dt2;
      node.um = u;
    }

    for (size_t i = 0 ; i < Fem::connected_elements.size() ; i++) {
      size_t ielem = Fem::connected_elements[i];
      Element& element = Fem::elements[ielem];
      Eigen::VectorXd u = Eigen::VectorXd::Zero(element.dof);

      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        u += element.nodes_p[inode]->u;
      }
      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        element.nodes_p[inode]->u = u/element.nnode;
      }
    }

    for (size_t i = 0 ; i < Fem::output_elements.size() ; i++) {
      size_t ielem = Fem::output_elements[i];
      Element& element = Fem::elements[ielem];
      element.calc_stress();
    }
  }

// ------------------------------------------------------------------- //
