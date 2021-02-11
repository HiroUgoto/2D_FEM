class Fem {
public:
  size_t nnode, nelem, dof;
  std::vector<Node> nodes;
  std::vector<Element> elements;
  std::vector<Material> materials;

  std::vector<size_t> free_nodes, fixed_nodes;
  std::vector<size_t> input_elements;
  std::vector<size_t> connected_elements;

  Fem (size_t dof, std::vector<Node> nodes,
                std::vector<Element> elements,
                std::vector<Material> materials);

  void set_init();

private:
  void _set_mesh();
  void _set_initial_condition();
  void _set_initial_matrix();
};
