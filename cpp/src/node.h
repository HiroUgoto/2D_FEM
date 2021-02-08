class Node {
public:
  int id;
  int dof;
  std::vector<double> xyz;
  std::vector<int> freedom;

  Node (int id, std::vector<double> xyz, std::vector<int> freedom);
  void print() ;
};
