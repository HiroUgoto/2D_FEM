class Element {
public:
  size_t id, material_id;
  std::string style;
  std::vector<size_t> inode;
  double gravity;

  std::vector<Node*> nodes_p;
  std::vector<Eigen::VectorXd*> u_p, v_p;

  Material* material_p;
  double rho, mass;

  size_t nnode, dim, ng, ng_all, dof, ndof;
  Eigen::VectorXd xi, w;
  Eigen::MatrixXd xnT;
  std::vector<Eigen::VectorXd> n_list;
  std::vector<Eigen::MatrixXd> dn_list;
  std::vector<double> w_list;

  Eigen::VectorXd M_diag, K_diag, C_diag;
  Eigen::MatrixXd K, K_off_diag, C, C_off_diag;
  Eigen::MatrixXd De, Dv;
  Eigen::VectorXd force;


  Element (size_t id, std::string style, size_t material_id, std::vector<size_t> inode);
  void print() ;

private:
  void set_style();

public:
  void set_nodes(std::vector<Node*> nodes_p);
  void set_material(Material* material_p);
  void set_pointer_list();
  void set_xn();

  void mk_local_matrix_init(const size_t dof);
  void mk_local_matrix();
  void mk_local_vector();


private:
  Eigen::MatrixXd mk_m(const Eigen::MatrixXd N);
  Eigen::MatrixXd mk_n(const size_t dof, const size_t nnode, const Eigen::VectorXd n);

  Eigen::MatrixXd mk_k(const Eigen::MatrixXd B, const Eigen::MatrixXd D);
  Eigen::MatrixXd mk_b(const size_t dof, const size_t nnode, const Eigen::MatrixXd dnj);
  Eigen::MatrixXd mk_b_T(const size_t dof, const size_t nnode, const Eigen::MatrixXd dnj);

  std::tuple<double, Eigen::MatrixXd> mk_dnj(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn);
  std::tuple<double, Eigen::Matrix2d> mk_inv_jacobi(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn);
  std::tuple<double, Eigen::Matrix2d> mk_jacobi(const Eigen::MatrixXd xnT, const Eigen::MatrixXd dn);


};
