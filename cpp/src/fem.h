class Fem {
public:
  int dof;
  std::vector<Node> nodes;
  std::vector<Element> elements;
  std::vector<Material> materials;

  Fem (int dof, std::vector<Node> nodes,
                std::vector<Element> elements,
                std::vector<Material> materials);
};
