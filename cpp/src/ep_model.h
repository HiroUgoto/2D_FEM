#include <unsupported/Eigen/CXX11/Tensor>

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;

class StateParameters {
  public:
    EM strain, stress;
    EM dstress, dstrain;

    double pmin;
    double p, R, dp;
    EM sij, rij;

    double Ge ,Ke;
    double psi, g;
    double Kp1, Kp1_b, Kp2, Kp2_b, D1, D2, h;
    double R_bar, g_bar, rho1_ratio, rho2_ratio, p_bar;
    EM rij_bar;

    EM nij, mij;
    double B;
    EM Tij, Zij;

    bool elastic_flag1, elastic_flag2;

    StateParameters(EM strain, EM stress, EM dstrain, EM dstress, bool ef1, bool ef2);

    void set_stress_variable();
    void set_stress_increment();
};


class Li2002: public EP {
  public:
    double G0,nu;
    double M,c,eg,rlambdac,xi;
    double d1,m,h1,h2,h3,n;
    double d2,h4;
    double a,b1,b2,b3;

    double eps;
    double pr,pmin;

    EM alpha;
    double beta,H1,H2;
    double L1;

    EM Z3,I3;

    double sqrt2_3,fn,rlambda_coeff,G2_coeff,g0,dg0;

    Li2002 (double G0=125, double nu=0.25, double M=1.25, double eg=0.934, double d1=0.41,
            double c=0.75, double rlambdac=0.019, double xi=0.7,
            double m=3.5, double h1=3.15, double h2=3.05, double h3=2.2, double n=1.1,
            double d2=1, double h4=3.5, double a=1, double b1=0.005, double b2=2, double b3=0.01);

    void print();
    void clear_strain();

    void initial_state(EV init_stress);
    std::tuple<EM, EV> set_Dp_matrix(EV FEMdstrain);
    EV strain_to_stress(EV FEMdstrain);

    void isotropic_compression(const double e, const double compression_stress);

    std::tuple<EM, StateParameters> plastic_deformation_strain(const EM dstrain_given, const EM dstress_given, const StateParameters sp0);

    std::tuple<EM, StateParameters> plastic_deformation_stress(const EM dstrain_given, const EM dstress_given, const StateParameters sp0);

  private:
    std::tuple<double, double> set_stress_variable(const EM stress);
    double set_stress_variable_p(const EM stress);
    std::tuple<double, double> set_strain_variable(const EM strain);
    double set_strain_variable_ev(const EM strain);

    Eigen::Tensor<double,4> isotropic_compression_stiffness(const double e, const double p);
    Eigen::Tensor<double,4> elastic_stiffness(const double G);
    double elastic_modulus_G(const double e, const double p);

  public:
    std::tuple<double, double> elastic_modulus(const double e, const double p);
    std::tuple<double, double> elastic_modulus_lame();

  private:
    EM solve_strain(const EM stress_mat, const Eigen::Tensor<double,4> E);
    EM solve_stress(const EM strain_mat, const Eigen::Tensor<double,4> E);

    std::tuple<double, double> check_unload(StateParameters sp);
    Eigen::Tensor<double,4> plastic_stiffness(StateParameters &sp);
    void update_parameters(StateParameters &sp);

    void set_mapping_stress(StateParameters &sp);
    std::tuple<double, EM, double, double> _F1_boundary_surface_all(const double t, const EM rij, const EM alpha);
    double _F1_boundary_surface(const double t, const EM rij, const EM alpha);
    std::tuple<EM, double, double> _mapping_r(const double t, const EM rij, const EM alpha);
    std::tuple<double, double> _mapping_r_Rg(const double t, const EM rij, const EM alpha);
    double _find_rho1_ratio(const EM rij, const EM alpha);


    void set_parameters(StateParameters &sp);
    double _accumulated_load_index(const double L1);
    double _scaling_factor(const double e, const double rho1_ratio);
    std::tuple<double, double> _plastic_modulus1(const double G, const double R_bar, const double g_bar, const double rho1_ratio, const double h, const double psi);
    std::tuple<double, double> _plastic_modulus2(const double G, const double Mg_R, const double rho2_ratio, const double sign);
    double _dilatancy1(const double R, const double g, const double rho1_ratio, const double psi);
    double _dilatancy2(const double Mg_R, const double sign);

    void set_parameter_nm(StateParameters &sp);
    EM _dF1_r(const EM rij_bar, const double R_bar, const double theta3_bar, const double g_bar, const double dg_bar);

    void set_parameter_TZ(StateParameters &sp);

    Eigen::Tensor<double,4> set_tensor_Ep(StateParameters &sp);

    double state_parameter(const double e, const double p);
    double Lode_angle(const EM dev_stress);
    std::tuple<double, double> Lode_angle_J2(const EM dev_stress);
    double g_theta(const EM dev_stress);
    std::tuple<double, double> g_theta_J2(const EM dev_stress);
    double dg_theta(const double theta3);

    Eigen::Tensor<double,4> set_Dijkl();
    Eigen::Tensor<double,4> set_Dikjl();
    Eigen::Tensor<double,4> set_Diljk();
};
