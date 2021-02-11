class Material {
public:
  size_t id;
  std::string style;
  std::vector<double> param;
  double rmu, rlambda, rho;

  Material (size_t id, std::string style, std::vector<double> param);
  void print() ;

private:
  void set_param();

public:
  Eigen::MatrixXd mk_d(const size_t dof);
  Eigen::MatrixXd mk_visco(const size_t dof);
  Eigen::MatrixXd mk_imp(const size_t dof);

};
