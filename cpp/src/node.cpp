#include "all.h"
#include <Eigen/Core>
#include "node.h"


Node::Node () {}
Node::Node (size_t id, std::vector<double> xyz, std::vector<size_t> freedom)
  : xyz(2)
  {
    Node::id = id;
    Node::xyz = xyz;
    Node::freedom = freedom;
    Node::dof = freedom.size();
  }

void Node::print() {
  std::cout << Node::id << ": ";
  for (size_t i = 0 ; i < Node::xyz.size() ; i++) {
    std::cout << Node::xyz.at(i) << " ";
  }
  std::cout << "-- ";
  for (size_t i = 0 ; i < Node::freedom.size() ; i++) {
    std::cout << Node::freedom.at(i) << " ";
  }
  std::cout << "\n";
}

void Node::set_initial_condition() {
  Node::u  = Eigen::VectorXd::Zero(Node::dof);
  Node::um = Eigen::VectorXd::Zero(Node::dof);
  Node::v  = Eigen::VectorXd::Zero(Node::dof);

  Node::mass = Eigen::VectorXd::Zero(Node::dof);
  Node::c    = Eigen::VectorXd::Zero(Node::dof);
  Node::k    = Eigen::VectorXd::Zero(Node::dof);

  Node::force = Eigen::VectorXd::Zero(Node::dof);
  Node::static_force  = Eigen::VectorXd::Zero(Node::dof);
  Node::dynamic_force = Eigen::VectorXd::Zero(Node::dof);

  Node::inv_mc = Eigen::VectorXd::Zero(Node::dof);
  Node::mass_inv_mc = Eigen::VectorXd::Zero(Node::dof);
  Node::c_inv_mc = Eigen::VectorXd::Zero(Node::dof);
  Node::dtdt_inv_mc = Eigen::VectorXd::Zero(Node::dof);

  Node::_up = Eigen::VectorXd::Zero(Node::dof);
  Node::_ur = Eigen::VectorXd::Zero(Node::dof);
  Node::_uy = Eigen::VectorXd::Zero(Node::dof);

}
