using EV = Eigen::VectorXd ;

class Node {
  public:
    size_t id, dof;
    std::vector<double> xyz;
    std::vector<size_t> freedom;

    EV u, um, v, u0;
    EV mass, c, k;
    EV force, static_force, dynamic_force;

    EV inv_mc, mass_inv_mc, c_inv_mc, dtdt_inv_mc;
    EV _up, _ur, _uy;


    Node ();
    Node (size_t id, std::vector<double> xyz, std::vector<size_t> freedom);

    void print();
    void set_initial_condition();
};
