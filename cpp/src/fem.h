using EV = Eigen::VectorXd ;

class Fem {
  public:
    size_t nnode, nelem, dof;
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<Material> materials;

    std::vector<Node*> free_nodes_p, fixed_nodes_p;
    std::vector<Element*> connected_elements_p;
    std::vector<Element*> input_elements_p;

    size_t output_nnode, output_nelem;
    std::vector<Node*> output_nodes_p;
    std::vector<Element*> output_elements_p;

    double dt, inv_dt2;
    

    Fem (size_t dof, std::vector<Node> nodes,
                  std::vector<Element> elements,
                  std::vector<Material> materials);

  public:
    void set_init();

  private:
    void _set_mesh();
    void _set_initial_condition();
    void _set_initial_matrix();

  public:
    void set_output(std::tuple<std::vector<size_t>, std::vector<size_t>> outputs);

  public:
    void self_gravity();

  private:
    void _self_gravity_cg(const bool full=true);

  public:
    void update_init(const double dt);

    void update_time(const EV acc0, const EV vel0, const bool input_wave=false, const bool FD=false);
    void update_matrix();
};
