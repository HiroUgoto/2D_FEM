class Element {
public:
  size_t id, material_id;
  std::string style;
  std::vector<size_t> inode;
  double gravity;

  size_t nnode, dim, ng;
  ElementStyle* estyle_p;
  Eigen::VectorXd xi, w;

  Element (size_t id, std::string style, size_t material_id, std::vector<size_t> inode);
  void print() ;

private:
  void set_style();
};
