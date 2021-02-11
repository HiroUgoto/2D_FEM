class Node {
  public:
    size_t id, dof;
    std::vector<double> xyz;
    std::vector<size_t> freedom;

    Eigen::VectorXd u, um, v, u0;
    Eigen::VectorXd mass, c, k;
    Eigen::VectorXd force, static_force, dynamic_force;

    Eigen::VectorXd inv_mc, mass_inv_mc, c_inv_mc, dtdt_inv_mc;
    Eigen::VectorXd _up, _ur, _uy;


    Node ();
    Node (size_t id, std::vector<double> xyz, std::vector<size_t> freedom);

    void
      print();
    void
      set_initial_condition();
};
