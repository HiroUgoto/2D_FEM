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
    for (auto& element : this->elements) {

      std::vector<Node*> nodes_p(element.nnode);
      for (size_t i = 0 ; i < element.nnode ; i++ ){
        size_t inode = element.inode[i];

        if (this->nodes[inode].id == inode) {
          nodes_p[i] = &this->nodes[inode];
        } else {
          for (auto& node : this->nodes) {
            if (node.id == inode) {
              nodes_p[i] = &node;
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
        for (auto& material : this->materials) {
          if (material.id == element.material_id) {
            material_p = &material;
            break;
          }
        }
      }

      element.set_material(material_p);

      if (element.style.find("input") != std::string::npos) {
        this->input_elements_p.push_back(&element);
      }
      if (element.style.find("connect") != std::string::npos) {
        this->connected_elements_p.push_back(&element);
      }
    }
  }

// ----------------------------- //
void Fem::_set_initial_condition() {
    for (auto& node : this->nodes) {
      node.set_initial_condition();

      size_t total_dof = std::accumulate(node.freedom.begin(),node.freedom.end(),0);
      if (total_dof == dof) {
        this->free_nodes_p.push_back(&node);
      } else {
        this->fixed_nodes_p.push_back(&node);
      }
    }
  }

// ----------------------------- //
void Fem::_set_initial_matrix(){
    for (auto& element : this->elements) {
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
    for (size_t inode = 0 ; inode < this->output_nnode ; inode++) {
      size_t id = output_node_list[inode];
      this->output_nodes_p.push_back(&this->nodes[id]);
    }

    this->output_nelem = output_element_list.size();
    for (size_t ielem = 0 ; ielem < this->output_nelem ; ielem++) {
      size_t id = output_element_list[ielem];
      this->output_elements_p.push_back(&this->elements[id]);
    }
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::self_gravity() {
    // Initial condition //
    double H = 0.0;
    for (auto& node : this->nodes) {
      if (node.xyz[1] > H) {
        H = node.xyz[1];
      }
    }

    double g = 9.8;
    double vp = 1500.0;
    for (auto& node : this->nodes) {
      node.u[0] = 0.0;
      node.u[1] = g/(2*vp*vp) * (H*H - node.xyz[1]*node.xyz[1]);
      node.um = node.u;
    }

    this->_self_gravity_cg(false);
    this->_self_gravity_cg(true);

    for (auto& node : this->nodes) {
      node.u0 = node.u;
      node.um = node.u;
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
    for (auto& node : this->nodes) {
      node.force = EV::Zero(node.dof);
    }
    for (auto& element : this->elements) {
      element.mk_ku();
    }
    for (auto& node : this->nodes) {
      for (size_t i = id ; i < node.dof ; i++) {
        if (node.freedom[i] == 0) {
          node._ur[i] = 0.0;
        } else {
          node._ur[i] = node.static_force[i] - node.force[i];
        }
      }
    }
    for (auto& element_p : this->connected_elements_p) {
      EV u = EV::Zero(element_p->dof);
      for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
        u += element_p->nodes_p[inode]->_ur;
      }
      for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
        element_p->nodes_p[inode]->_ur = u/element_p->nnode;
      }
    }
    for (auto& element_p : this->input_elements_p) {
      for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
        element_p->nodes_p[inode]->_ur = EV::Zero(element_p->dof);
      }
    }
    for (auto& node : this->nodes) {
      node._up = node._ur;
    }


    // --- CG iterations --- //
    for (size_t it=0 ; it < 10*this->nnode ; it++) {
      double rr, rr1, py, alpha, beta;

      // y = Ap
      for (auto& node : this->nodes) {
        node.force = EV::Zero(node.dof);
      }
      for (auto& element : this->elements) {
        element.mk_ku_up();
      }
      for (auto& node : this->nodes) {
        node._uy = node.force;
      }

      // correction boundary condition
      for (auto& node : this->nodes) {
        for (size_t i = id ; i < node.dof ; i++) {
          if (node.freedom[i] == 0) {
            node._uy[i] = 0.0;
          }
        }
      }
      for (auto& element_p : this->connected_elements_p) {
        EV u = EV::Zero(element_p->dof);
        for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
          u += element_p->nodes_p[inode]->_uy;
        }
        for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
          element_p->nodes_p[inode]->_uy = u/element_p->nnode;
        }
      }
      for (auto& element_p : this->input_elements_p) {
        for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
          element_p->nodes_p[inode]->_uy = EV::Zero(element_p->dof);
        }
      }

      // alpha = rr/py
      rr = 0.0; py = 0.0;
      for (auto& node : this->nodes) {
        rr += node._ur.dot(node._ur);
        py += node._up.dot(node._uy);
      }
      alpha = rr/py;

      // x = x + alpha*p
      rr1 = 0.0;
      for (auto& node : this->nodes) {
        for (size_t i = id ; i < node.dof ; i++) {
          if (node.freedom[i] != 0) {
            node.u[i] += alpha * node._up[i];
            node._ur[i] -= alpha * node._uy[i];
          }
        }
        rr1 += node._ur.dot(node._ur);
      }

      if (rr1 < 1.e-10) {
        std::cout << " (self gravity process .. ) " << it << " " ;
        std::cout << this->nodes[0].u[1] << " ";
        std::cout << rr1 << "\n";
        break;
      }

      // p = r + beta*p
      beta = rr1/rr;
      for (auto& node : this->nodes) {
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
    for (auto& node : this->nodes) {
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
    if(FD) {
      if(input_wave) {
        this->update_time_input_FD(vel0);
      } else {
        this->update_time_FD(acc0);
      }
    } else {
      if(input_wave) {
        this->update_time_input_MD(vel0);
      } else {
        this->update_time_MD(acc0);
      }
    }
  }

// ------------------------------------------------------------------- //
void Fem::update_time_FD(const EV acc0) {

    this->_update_matrix();

    for (auto& node : this->nodes) {
      node.force = -node.dynamic_force;
    }

    for (auto& element : this->elements) {
      element.update_bodyforce(acc0);
    }

    for (auto& element : this->elements) {
      element.mk_B_stress();
      element.mk_cv();
    }

    this->_update_time_set_free_nodes();
    this->_update_time_set_fixed_nodes();
    this->_update_time_set_connected_elements();

    for (auto& element_p : this->output_elements_p) {
      element_p->calc_stress();
    }
  }

// ------------------------------------------------------------------- //
void Fem::update_time_input_FD(const EV vel0) {

    this->_update_matrix();

    for (auto& node : this->nodes) {
      node.force = -node.dynamic_force;
    }

    for (auto& element_p : this->input_elements_p) {
      element_p->update_inputwave(vel0);
    }

    for (auto& element : this->elements) {
      element.mk_B_stress();
      element.mk_cv();
    }

    this->_update_time_set_free_nodes();
    this->_update_time_set_fixed_nodes();
    this->_update_time_set_connected_elements();

    for (auto& element_p : this->output_elements_p) {
      element_p->calc_stress();
    }
  }


// ------------------------------------------------------------------- //
void Fem::update_time_MD(const EV acc0) {
    for (auto& node : this->nodes) {
      node.force = -node.dynamic_force;
    }

    for (auto& element : this->elements) {
      element.update_bodyforce(acc0);
    }

    for (auto& element : this->elements) {
      element.mk_ku_cv();
    }

    this->_update_time_set_free_nodes();
    this->_update_time_set_fixed_nodes();
    this->_update_time_set_connected_elements();

    for (auto& element_p : this->output_elements_p) {
      element_p->calc_stress();
    }
  }

// ------------------------------------------------------------------- //
void Fem::update_time_input_MD(const EV vel0) {
    for (auto& node : this->nodes) {
      node.force = -node.dynamic_force;
    }

    for (auto& element_p : this->input_elements_p) {
      element_p->update_inputwave(vel0);
    }

    for (auto& element : this->elements) {
      element.mk_ku_cv();
    }

    this->_update_time_set_free_nodes();
    this->_update_time_set_fixed_nodes();
    this->_update_time_set_connected_elements();

    for (auto& element_p : this->output_elements_p) {
      element_p->calc_stress();
    }
  }

// ------------------------------------------------------------------- //
// ------------------------------------------------------------------- //
void Fem::_update_time_set_free_nodes() {
    for (auto& node_p : this->free_nodes_p) {
      EV u = node_p->u;
      for (size_t i = 0 ; i < node_p->dof ; i++) {
        node_p->u[i] = node_p->mass_inv_mc[i]*(2.0*u[i] - node_p->um[i])
                + node_p->c_inv_mc[i]*node_p->um[i] - node_p->dtdt_inv_mc[i]*node_p->force[i];
      }
      node_p->v = (node_p->u - node_p->um) * this->inv_dt2;
      node_p->um = u;
    }
  }

void Fem::_update_time_set_fixed_nodes() {
    for (auto& node_p : this->fixed_nodes_p) {
      EV u = node_p->u;
      for (size_t i = 0 ; i < node_p->dof ; i++) {
        if (node_p->freedom[i] == 0) {
          node_p->u[i] = 0.0;
        } else {
          node_p->u[i] = node_p->mass_inv_mc[i]*(2.0*u[i] - node_p->um[i])
                  + node_p->c_inv_mc[i]*node_p->um[i] - node_p->dtdt_inv_mc[i]*node_p->force[i];
        }
      }
      node_p->v = (node_p->u - node_p->um) * this->inv_dt2;
      node_p->um = u;
    }
  }

void Fem::_update_time_set_connected_elements() {
    for (auto& element_p : this->connected_elements_p) {
      EV u = EV::Zero(element_p->dof);
      for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
        u += element_p->nodes_p[inode]->u;
      }
      for (size_t inode = 0 ; inode < element_p->nnode ; inode++) {
        element_p->nodes_p[inode]->u = u/element_p->nnode;
      }
    }
  }

// ------------------------------------------------------------------- //
void Fem::_update_matrix() {
    for (auto& node : this->nodes) {
      node.mass = EV::Zero(node.dof);
      node.c    = EV::Zero(node.dof);
      node.dynamic_force = EV::Zero(node.dof);
    }

    for (auto& element : this-> elements) {
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

    for (auto& node : this->nodes) {
      for (size_t i = 0 ; i < node.dof ; i++) {
        node.inv_mc[i] = 1.0 / (node.mass[i] + 0.5*dt*node.c[i]);
        node.mass_inv_mc[i] = node.mass[i] * node.inv_mc[i];
        node.c_inv_mc[i] = node.c[i] * node.inv_mc[i] * 0.5*dt;
        node.dtdt_inv_mc[i] = dt*dt*node.inv_mc[i];
      }
    }
  }
