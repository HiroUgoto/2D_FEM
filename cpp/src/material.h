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
};
