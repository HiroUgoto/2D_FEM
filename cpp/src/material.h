using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;

class Material {
  public:
    size_t id;
    std::string style;
    std::vector<double> param, param_ep;
    double rmu, rlambda, rho;
    double kv, kh;
    EM R;

    Material();
    Material(size_t id, std::string style, std::vector<double> param);

    void set_init(size_t id, std::string style, std::vector<double> param);
    void print();

  private:
    void set_param();

  public:
    EM mk_d(const size_t dof);
    EM mk_d_spring();
    EM mk_visco(const size_t dof);
    EM mk_imp(const size_t dof);
};
