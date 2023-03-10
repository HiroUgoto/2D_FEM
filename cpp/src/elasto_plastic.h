#include <unsupported/Eigen/CXX11/Tensor>

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;

class EP {
  public:
    size_t dof;
    std::string style;
    double e0;

    EM stress,strain;
    double e;
    double psi, fL, h, n_e;
    double out_H1, out_H2, out_L1, out_h1, out_h2, out_D1, out_D2, out_Kp1, out_Kp2, out_Ge, out_Ke;

    EP ();
    virtual ~EP ();


    EM FEMstress_to_matrix(EV FEMstress);
    EM FEMstrain_to_matrix(EV FEMstrain);

    EV matrix_to_FEMstress(EM stress);
    EM modulus_to_Dmatrix(Eigen::Tensor<double,4> E);

    EM vector_to_matrix(EV vec);
    EV matrix_to_vector(EM mat);

    virtual void print() = 0;
    virtual void clear_strain() = 0;

    virtual std::tuple<double, double> elastic_modulus(const double e, const double p) = 0;
    virtual std::tuple<double, double> elastic_modulus_lame() = 0;
    virtual void initial_state_isotropic(EV init_stress) = 0;
    virtual void initial_state(EV init_stress) = 0;
    virtual std::tuple<EM, EV, double> set_Dp_matrix(EV FEMdstrain) = 0;
    virtual EV strain_to_stress(EV FEMdstrain) = 0;
};

// ------------------------------------------------------------------- //
EP* set_ep_style(double dof, std::string style, std::vector<double> param_ep);
