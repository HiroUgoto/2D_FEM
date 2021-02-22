using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;

class ElementStyle {
  public:
    size_t dim;
    size_t ng;
    EV xi, w;

    size_t ng_all;
    std::vector<EV> n_list;
    std::vector<EM> dn_list;
    std::vector<double> w_list;

    EM dn_center;

    ElementStyle ();

    virtual EV shape_function_n (double xi, double zeta);
    virtual EM shape_function_dn (double xi, double zeta);
};

// ----------------------------------------------------- //
ElementStyle* set_element_style(const std::string style);
void set_gauss_points (const size_t n, EV& xi, EV& w);

// ----------------------------------------------------- //
class Solid_2d_4Node: public ElementStyle {
  public:
    Solid_2d_4Node ();
    EV shape_function_n (double xi, double zeta);
    EM shape_function_dn (double xi, double zeta);
};

// ----------------------------------------------------- //
class Solid_2d_8Node: public ElementStyle {
  public:
    Solid_2d_8Node ();
    EV shape_function_n (double xi, double zeta);
    EM shape_function_dn (double xi, double zeta);
};

// ----------------------------------------------------- //
class Solid_2d_9Node: public ElementStyle {
  public:
    Solid_2d_9Node ();
    EV shape_function_n (double xi, double zeta);
    EM shape_function_dn (double xi, double zeta);
};

// ----------------------------------------------------- //
class Line_1d_2Node: public ElementStyle {
  public:
    Line_1d_2Node ();
    EV shape_function_n (double xi, double zeta=0.0);
    EM shape_function_dn (double xi, double zeta=0.0);
};

// ----------------------------------------------------- //
class Line_1d_3Node: public ElementStyle {
  public:
    Line_1d_3Node ();
    EV shape_function_n (double xi, double zeta=0.0);
    EM shape_function_dn (double xi, double zeta=0.0);
};

// ----------------------------------------------------- //
class Input_1d_2Node: public Line_1d_2Node {};
class Input_1d_3Node: public Line_1d_3Node {};

// ----------------------------------------------------- //
class Connect: public Line_1d_2Node {
  public:
    Connect ();
};
