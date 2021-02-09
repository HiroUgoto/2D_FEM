class Fem {
public:
  size_t dof;
  std::vector<Node> nodes;
  std::vector<Element> elements;
  std::vector<Material> materials;

  Fem (size_t dof, std::vector<Node> nodes,
                std::vector<Element> elements,
                std::vector<Material> materials);
};
