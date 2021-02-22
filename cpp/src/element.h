class Element {
  public:
    size_t id;
    int material_id;
    std::string style;
    std::vector<size_t> inode;
    double gravity;

    std::vector<Node*> nodes_p;
    std::vector<Eigen::VectorXd*> u_p, v_p;

    Material material;
    double rho, mass;

    size_t nnode, dim, ng, ng_all, dof, ndof;
    Eigen::VectorXd xi, w;
    Eigen::MatrixXd xnT;
    std::vector<Eigen::VectorXd> n_list;
    std::vector<Eigen::MatrixXd> dn_list;
    std::vector<double> w_list;
    Eigen::MatrixXd dn_center;

    Eigen::VectorXd M_diag, K_diag, C_diag;
    Eigen::MatrixXd K, K_off_diag, C, C_off_diag;
    Eigen::MatrixXd De, Dv, imp;
    Eigen::VectorXd force;

    Eigen::VectorXd strain, stress;

    std::vector<Eigen::VectorXd*> _up_p;

    Element (size_t id, std::string style, int material_id, std::vector<size_t> inode);

    void
      print() ;

  private:
    void
      set_style();

  public:
    void
      set_nodes(std::vector<Node*> nodes_p);
    void
      set_material(Material* material_p);
    void
      set_pointer_list();
    void
      set_xn();

    void
      mk_local_matrix_init(const size_t dof);
    void
      mk_local_matrix();
    void
      mk_local_vector();
    void
      mk_local_update();

    void
      mk_ku();
    void
      mk_ku_u(const std::vector<Eigen::VectorXd*> u_p);

    void
      mk_cv();
    void
      mk_ku_cv();

    void
      mk_B_stress();


  private:
    Eigen::VectorXd
      mk_u_hstack(const std::vector<Eigen::VectorXd*> u_p);
    Eigen::MatrixXd
      mk_u_vstack(const std::vector<Eigen::VectorXd*> u_p);

  public:
    void
      update_bodyforce(const Eigen::VectorXd acc0);
    void
      mk_bodyforce(const Eigen::VectorXd acc0);

    void
      update_inputwave(const Eigen::VectorXd vel0);


    void
      calc_stress();
};

Eigen::MatrixXd
  mk_m(const Eigen::MatrixXd N);
Eigen::MatrixXd
  mk_n(const size_t dof, const size_t nnode, const Eigen::VectorXd n);

Eigen::MatrixXd
  mk_nqn(const Eigen::MatrixXd N, const Eigen::MatrixXd q, const Eigen::MatrixXd imp);
std::tuple<double, Eigen::MatrixXd>
  mk_q(const size_t dof, const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn);

Eigen::MatrixXd
  mk_k(const Eigen::MatrixXd B, const Eigen::MatrixXd D);
Eigen::MatrixXd
  mk_b(const size_t dof, const size_t nnode, const Eigen::MatrixXd dnj);
Eigen::MatrixXd
  mk_b_T(const size_t dof, const size_t nnode, const Eigen::MatrixXd dnj);

Eigen::VectorXd
  Hencky_stress(const Eigen::MatrixXd D, const Eigen::MatrixXd dnj, const Eigen::MatrixXd u);
std::tuple<double, Eigen::Matrix2d>
  Euler_log_strain(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u);
std::tuple<double, Eigen::Matrix2d>
  mk_FF(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u);
std::tuple<double, Eigen::Matrix2d>
  mk_F(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u);
Eigen::Matrix2d
  mk_dnu(const Eigen::MatrixXd dnj, const Eigen::MatrixXd u);

std::tuple<double, Eigen::MatrixXd>
  mk_dnj(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn);
std::tuple<double, Eigen::Matrix2d>
  mk_inv_jacobi(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn);
std::tuple<double, Eigen::Matrix2d>
  mk_jacobi(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn);
