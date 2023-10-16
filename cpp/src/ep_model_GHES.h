#include <unsupported/Eigen/CXX11/Tensor>

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;

class GHE_S: public EP {
    public:
        double G0,nu;
        double gr,hmax;

        GHE_S (double G0, double gr, double hmax, double nu);

        void print();
        void clear_strain();

        std::tuple<double, double> elastic_modulus(const double e, const double p);
        std::tuple<double, double> elastic_modulus_lame();

        void initial_state_isotropic(EV init_stress);
        void initial_state(EV init_stress);

        std::tuple<EM, EV, double> set_Dp_matrix(EV FEMdstrain);
        EV strain_to_stress(EV FEMdstrain);


};