#include "all.h"
#include "node.h"
#include <Eigen/Core>


Node::Node (int id, std::vector<double> xyz, std::vector<int> freedom)
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
