#include "all.h"
#include <Eigen/Core>
#include "node.h"

using EV = Eigen::VectorXd ;

Node::Node () {}
Node::Node (size_t id, std::vector<double> xyz, std::vector<size_t> freedom)
  {
    this->id = id;
    this->xyz = xyz;
    this->freedom = freedom;
    this->dof = freedom.size();
  }

void Node::print() {
  std::cout << this->id << ": ";
  for (size_t i = 0 ; i < this->xyz.size() ; i++) {
    std::cout << this->xyz.at(i) << " ";
  }
  std::cout << "-- ";
  for (size_t i = 0 ; i < this->freedom.size() ; i++) {
    std::cout << this->freedom.at(i) << " ";
  }
  std::cout << "\n";
}

void Node::set_initial_condition() {
  this->u = EV::Zero(this->dof);
  this->um = EV::Zero(this->dof);
  this->v  = EV::Zero(this->dof);

  this->mass = EV::Zero(this->dof);
  this->c    = EV::Zero(this->dof);
  this->k    = EV::Zero(this->dof);

  this->force = EV::Zero(this->dof);
  this->static_force  = EV::Zero(this->dof);
  this->dynamic_force = EV::Zero(this->dof);

  this->inv_mc = EV::Zero(this->dof);
  this->mass_inv_mc = EV::Zero(this->dof);
  this->c_inv_mc = EV::Zero(this->dof);
  this->dtdt_inv_mc = EV::Zero(this->dof);

  this->_up = EV::Zero(this->dof);
  this->_ur = EV::Zero(this->dof);
  this->_uy = EV::Zero(this->dof);

}
