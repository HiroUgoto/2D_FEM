class Node {
public:
  size_t id, dof;
  std::vector<double> xyz;
  std::vector<size_t> freedom;

  Node (size_t id, std::vector<double> xyz, std::vector<size_t> freedom);
  void print() ;
};
