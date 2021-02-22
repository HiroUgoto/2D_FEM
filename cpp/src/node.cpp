#include "all.h"
#include <Eigen/Core>
#include "node.h"


Node::Node () {}
Node::Node (size_t id, std::vector<double> xyz, std::vector<size_t> freedom)
  : xyz(2)
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
  this->u  = Eigen::VectorXd::Zero(this->dof);
  this->um = Eigen::VectorXd::Zero(this->dof);
  this->v  = Eigen::VectorXd::Zero(this->dof);

  this->mass = Eigen::VectorXd::Zero(this->dof);
  this->c    = Eigen::VectorXd::Zero(this->dof);
  this->k    = Eigen::VectorXd::Zero(this->dof);

  this->force = Eigen::VectorXd::Zero(this->dof);
  this->static_force  = Eigen::VectorXd::Zero(this->dof);
  this->dynamic_force = Eigen::VectorXd::Zero(this->dof);

  this->inv_mc = Eigen::VectorXd::Zero(this->dof);
  this->mass_inv_mc = Eigen::VectorXd::Zero(this->dof);
  this->c_inv_mc = Eigen::VectorXd::Zero(this->dof);
  this->dtdt_inv_mc = Eigen::VectorXd::Zero(this->dof);

  this->_up = Eigen::VectorXd::Zero(this->dof);
  this->_ur = Eigen::VectorXd::Zero(this->dof);
  this->_uy = Eigen::VectorXd::Zero(this->dof);

}
