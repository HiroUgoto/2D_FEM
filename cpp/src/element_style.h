class ElementStyle {
public:
  size_t dim;
  size_t ng;
  Eigen::VectorXd xi, w;

  ElementStyle ();
  virtual Eigen::VectorXd shape_function_n (double xi, double zeta);
};

// ----------------------------------------------------- //
ElementStyle* set_element_style(const std::string style);

void set_gauss_points (const size_t n, Eigen::VectorXd& xi, Eigen::VectorXd& w);

// ----------------------------------------------------- //
class Solid_2d_4Node: public ElementStyle {
public:
  Solid_2d_4Node ();
  Eigen::VectorXd shape_function_n (double xi, double zeta) override;
};

// ----------------------------------------------------- //
class Solid_2d_9Node: public ElementStyle {
public:
  Solid_2d_9Node ();
  Eigen::VectorXd shape_function_n (double xi, double zeta) override;
};
