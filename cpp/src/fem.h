class Fem {
  public:
    size_t nnode, nelem, dof;
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<Material> materials;

    std::vector<size_t> free_nodes, fixed_nodes;
    std::vector<size_t> input_elements;
    std::vector<size_t> connected_elements;

    size_t output_nnode, output_nelem;
    std::vector<size_t> output_nodes, output_elements;

    double dt, inv_dt2;


    Fem (size_t dof, std::vector<Node> nodes,
                  std::vector<Element> elements,
                  std::vector<Material> materials);

  public:
    void
      set_init();

  private:
    void
      _set_mesh();
    void
      _set_initial_condition();
    void
      _set_initial_matrix();

  public:
    void
      set_output(std::tuple<std::vector<size_t>, std::vector<size_t>> outputs);

    void
      self_gravity();

  private:
    void
      _self_gravity_cg(const bool full=true);

  public:
    void
      update_init(const double dt);

    void
      update_time(const Eigen::VectorXd acc0, const Eigen::VectorXd vel0, const bool input_wave=false);

};
