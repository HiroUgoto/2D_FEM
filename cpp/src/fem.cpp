#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "fem.h"

using EV = Eigen::VectorXd ;

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
Fem::Fem (size_t dof, std::vector<Node> nodes,
                std::vector<Element> elements,
                std::vector<Material> materials)
  {
    this->nnode = nodes.size();
    this->nelem = elements.size();
    this->dof = dof;

    this->nodes = nodes;
    this->elements = elements;
    this->materials = materials;
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::set_init() {
    this->_set_mesh();
    this->_set_initial_condition();
    this->_set_initial_matrix();
  }

// ----------------------------- //
void Fem::_set_mesh() {
    for (size_t ielem = 0 ; ielem < this->nelem ; ielem++ ){
      Element& element = this->elements[ielem];

      std::vector<Node*> nodes_p(element.nnode);
      for (size_t i = 0 ; i < element.nnode ; i++ ){
        size_t inode = element.inode[i];

        if (this->nodes[inode].id == inode) {
          nodes_p[i] = &this->nodes[inode];
        } else {
          for (size_t i0 = 0 ; i0 < this->nnode ; i0++ ) {
            if (this->nodes[i0].id == inode) {
              nodes_p[i] = &this->nodes[i0];
              break;
            }
          }
        }
      }

      element.set_nodes(nodes_p);

      Material* material_p = nullptr;
      if (this->materials[element.material_id].id == element.material_id) {
        material_p = &this->materials[element.material_id];

      } else {
        for (size_t i0 = 0 ; i0 < this->materials.size() ; i0++) {
          if (this->materials[i0].id == element.material_id) {
            material_p = &this->materials[i0];
            break;
          }
        }
      }

      element.set_material(material_p);

      if (element.style.find("input") != std::string::npos) {
        this->input_elements.push_back(ielem);
      }
      if (element.style.find("connect") != std::string::npos) {
        this->connected_elements.push_back(ielem);
      }
    }
  }

// ----------------------------- //
void Fem::_set_initial_condition() {
    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];

      node.set_initial_condition();

      size_t total_dof = std::accumulate(node.freedom.begin(),node.freedom.end(),0);
      if (total_dof == dof) {
        this->free_nodes.push_back(inode);
      } else {
        this->fixed_nodes.push_back(inode);
      }
    }

    for (size_t ielem = 0 ; ielem < this->nelem ; ielem++ ){
      Element& element = this->elements[ielem];
      element.set_pointer_list();
    }
  }

// ----------------------------- //
void Fem::_set_initial_matrix(){
    for (size_t ielem = 0 ; ielem < this->nelem ; ielem++ ){
      Element& element = this->elements[ielem];

      element.set_xn();
      element.mk_local_matrix_init(this->dof);
      element.mk_local_matrix();
      element.mk_local_vector();

      size_t id = 0;
      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        for (size_t i = 0 ; i < this->dof ; i++) {
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
    auto [output_node_list, output_element_list] = outputs;

    this->output_nnode = output_node_list.size();
    this->output_nodes = output_node_list;

    this->output_nelem = output_element_list.size();
    this->output_elements = output_element_list;
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::self_gravity() {
    // Initial condition //
    double H = 0.0;
    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      if (this->nodes[inode].xyz[1] > H) {
        H = this->nodes[inode].xyz[1];
      }
    }

    double g = 9.8;
    double vp = 1500.0;
    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];

      node.u[0] = 0.0;
      node.u[1] = g/(2*vp*vp) * (H*H - node.xyz[1]*node.xyz[1]);
      node.um = node.u;
    }

    this->_self_gravity_cg(false);
    this->_self_gravity_cg(true);

    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];
      node.u0 = node.u;
    }

  }

// ------------------------------------------------------------------- //
void Fem::_self_gravity_cg(const bool full) {
    size_t id;
    if (full) {
      id = 0;
    } else {
      id = 1;
    }

    // CG Method //
    // --- set initial variables --- //
    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];
      node.force = EV::Zero(node.dof);
    }
    for (size_t ielem = 0 ; ielem < this->nelem ; ielem++) {
      Element& element = this->elements[ielem];
      element.mk_ku();
    }

    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];
      for (size_t i = id ; i < node.dof ; i++) {
        if (node.freedom[i] == 0) {
          node._ur[i] = 0.0;
        } else {
          node._ur[i] = node.static_force[i] - node.force[i];
        }
      }
    }
    for (size_t i = 0 ; i < this->connected_elements.size() ; i++) {
      size_t ielem = this->connected_elements[i];
      Element& element = this->elements[ielem];
      EV u = EV::Zero(element.dof);

      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        u += element.nodes_p[inode]->_ur;
      }
      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        element.nodes_p[inode]->_ur = u/element.nnode;
      }
    }

    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];
      node._up = node._ur;
    }


    // --- CG iterations --- //
    for (size_t it=0 ; it < 10*this->nnode ; it++) {
      double rr, rr1, py, alpha, beta;

      // y = Ap
      for (size_t inode = 0 ; inode < this->nnode ; inode++) {
        Node& node = this->nodes[inode];
        node.force = EV::Zero(node.dof);
      }
      for (size_t ielem = 0 ; ielem < this->nelem ; ielem++) {
        Element& element = this->elements[ielem];
        element.mk_ku_up();
      }
      for (size_t inode = 0 ; inode < this->nnode ; inode++) {
        Node& node = this->nodes[inode];
        node._uy = node.force;
      }

      // correction boundary condition
      for (size_t inode = 0 ; inode < this->nnode ; inode++) {
        Node& node = this->nodes[inode];
        for (size_t i = id ; i < node.dof ; i++) {
          if (node.freedom[i] == 0) {
            node._uy[i] = 0.0;
          }
        }
      }
      for (size_t i = 0 ; i < this->connected_elements.size() ; i++) {
        size_t ielem = this->connected_elements[i];
        Element& element = this->elements[ielem];
        EV u = EV::Zero(element.dof);

        for (size_t inode = 0 ; inode < element.nnode ; inode++) {
          u += element.nodes_p[inode]->_uy;
        }
        for (size_t inode = 0 ; inode < element.nnode ; inode++) {
          element.nodes_p[inode]->_uy = u/element.nnode;
        }
      }

      // alpha = rr/py
      rr = 0.0; py = 0.0;
      for (size_t inode = 0 ; inode < this->nnode ; inode++) {
        Node& node = this->nodes[inode];
        rr += node._ur.dot(node._ur);
        py += node._up.dot(node._uy);
      }
      alpha = rr/py;

      // x = x + alpha*p
      rr1 = 0.0;
      for (size_t inode = 0 ; inode < this->nnode ; inode++) {
        Node& node = this->nodes[inode];
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
      for (size_t inode = 0 ; inode < this->nnode ; inode++) {
        Node& node = this->nodes[inode];
        for (size_t i = id ; i < node.dof ; i++) {
          if (node.freedom[i] != 0) {
            node._up[i] = node._ur[i] + beta * node._up[i];
          }
        }
      }

      if (it%100 == 0){
        std::cout << " (self gravity process .. ) " << it << " " ;
        std::cout << this->nodes[0].u[1] << " ";
        std::cout << rr1 << "\n";
      }
    }

  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::update_init(const double dt) {
    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];

      for (size_t i = 0 ; i < node.dof ; i++) {
        node.inv_mc[i] = 1.0 / (node.mass[i] + 0.5*dt*node.c[i]);
        node.mass_inv_mc[i] = node.mass[i] * node.inv_mc[i];
        node.c_inv_mc[i] = node.c[i] * node.inv_mc[i] * 0.5*dt;
        node.dtdt_inv_mc[i] = dt*dt*node.inv_mc[i];
      }
    }

    this->dt = dt;
    this->inv_dt2 = 1.0/(2.0*dt);
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::update_time(const EV acc0, const EV vel0, const bool input_wave, const bool FD) {
    if (FD) {
      this->update_matrix();

    } else {
      for (size_t inode = 0 ; inode < this->nnode ; inode++) {
        Node& node = this->nodes[inode];
        node.dynamic_force = node.static_force;
      }
    }

    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];
      node.force = -node.dynamic_force;
    }

    if (input_wave) {
      for (size_t i = 0 ; i < this->input_elements.size() ; i++) {
        size_t ielem = this->input_elements[i];
        Element& element = this->elements[ielem];
        element.update_inputwave(vel0);
      }

    } else {
      for (size_t ielem = 0 ; ielem < this->nelem ; ielem++) {
        Element& element = this->elements[ielem];
        element.update_bodyforce(acc0);
      }
    }

    if (FD) {
      for (size_t ielem = 0 ; ielem < this->nelem ; ielem++) {
        Element& element = this->elements[ielem];
        element.mk_B_stress();
        element.mk_cv();
      }

    } else {
      for (size_t ielem = 0 ; ielem < this->nelem ; ielem++) {
        Element& element = this->elements[ielem];
        element.mk_ku_cv();
      }
    }

    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];
      EV u = node.u;

      for (size_t i = 0 ; i < node.dof ; i++) {
        if (node.freedom[i] == 0) {
          node.u[i] = 0.0;
        } else {
          node.u[i] = node.mass_inv_mc[i]*(2.0*u[i]-node.um[i]) + node.c_inv_mc[i]*node.um[i] - node.dtdt_inv_mc[i]*node.force[i];
        }
      }
      node.v = (node.u - node.um) * this->inv_dt2;
      node.um = u;
    }

    for (size_t i = 0 ; i < this->connected_elements.size() ; i++) {
      size_t ielem = this->connected_elements[i];
      Element& element = this->elements[ielem];
      EV u = EV::Zero(element.dof);

      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        u += element.nodes_p[inode]->u;
      }
      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        element.nodes_p[inode]->u = u/element.nnode;
      }
    }

    for (size_t i = 0 ; i < this->output_elements.size() ; i++) {
      size_t ielem = this->output_elements[i];
      Element& element = this->elements[ielem];
      element.calc_stress();
    }
  }

// ------------------------------------------------------------------- //
void Fem::update_matrix() {
    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];

      node.mass = EV::Zero(node.dof);
      node.c    = EV::Zero(node.dof);
      node.dynamic_force = EV::Zero(node.dof);
    }

    for (size_t ielem = 0 ; ielem < this->nelem ; ielem++) {
      Element& element = this->elements[ielem];

      element.set_xn();
      element.mk_local_update();

      size_t id = 0;
      for (size_t inode = 0 ; inode < element.nnode ; inode++) {
        for (size_t i = 0 ; i < this->dof ; i++) {
          element.nodes_p[inode]->mass[i] += element.M_diag[id];
          element.nodes_p[inode]->c[i] += element.C_diag[id];
          element.nodes_p[inode]->dynamic_force[i] += element.force[id];
          id++;
        }
      }
    }

    for (size_t inode = 0 ; inode < this->nnode ; inode++) {
      Node& node = this->nodes[inode];

      for (size_t i = 0 ; i < node.dof ; i++) {
        node.inv_mc[i] = 1.0 / (node.mass[i] + 0.5*dt*node.c[i]);
        node.mass_inv_mc[i] = node.mass[i] * node.inv_mc[i];
        node.c_inv_mc[i] = node.c[i] * node.inv_mc[i] * 0.5*dt;
        node.dtdt_inv_mc[i] = dt*dt*node.inv_mc[i];
      }
    }

  }
