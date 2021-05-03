using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;
using EM2 = Eigen::Matrix2d ;

class Element {
  public:
    size_t id;
    int material_id;
    std::string style;
    std::vector<size_t> inode;
    double gravity;

    std::vector<Node*> nodes_p;

    Material material;
    double rho, mass;

    size_t nnode, dim, ng, ng_all, dof, ndof;
    EV xi, w;
    EM xnT;
    std::vector<EV> n_list;
    std::vector<EM> dn_list;
    std::vector<double> w_list;
    EM dn_center;

    EV M_diag, K_diag, C_diag;
    EM K, K_off_diag, C, C_off_diag;
    EM De, Dv, imp;
    EV force;
    EV strain, stress;

    Element (size_t id, std::string style, int material_id, std::vector<size_t> inode);
    void print();

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
    void mk_local_update();

    void mk_ku();
    void mk_ku_up();
    void mk_cv();
    void mk_ku_cv();

    void mk_B_stress();

  private:
    EV mk_u_hstack();
    EV mk_v_hstack();
    EV mk_u_hstack_up();
    EM mk_u_vstack();

  public:
    void update_bodyforce(const EV acc0);
    void mk_bodyforce(const EV acc0);
    void update_inputwave(const EV vel0);
    void calc_stress();
};

EM mk_m(const EM N);
EM mk_n(const size_t dof, const size_t nnode, const EV n);

EM mk_nqn(const EM N, const EM q, const EM imp);
std::tuple<double, EM> mk_q(const size_t dof, const EM xnT, const EM dn);

EM mk_k(const EM B, const EM D);
EM mk_b(const size_t dof, const size_t nnode, const EM dnj);
EM mk_b_T(const size_t dof, const size_t nnode, const EM dnj);

EV Hencky_stress(const EM D, const EM dnj, const EM u);
std::tuple<double, EM2> Euler_log_strain(const EM dnj, const EM u);
std::tuple<double, EM2> mk_FF(const EM dnj, const EM u);
std::tuple<double, EM2> mk_F_inv(const EM dnj, const EM u);
std::tuple<double, EM2> mk_F(const EM dnj, const EM u);
EM2 mk_dnu(const EM dnj, const EM u);

std::tuple<double, EM>  mk_dnj(const EM xnT, const EM dn);
std::tuple<double, EM2> mk_inv_jacobi(const EM xnT, const EM dn);
std::tuple<double, EM2> mk_jacobi(const EM xnT, const EM dn);
