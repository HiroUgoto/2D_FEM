#include "all.h"
#include <Eigen/Core>
#include "node.h"


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
