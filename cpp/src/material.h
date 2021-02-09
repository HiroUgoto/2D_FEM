class Material {
public:
  int id;
  std::string style;
  std::vector<double> param;
  double rmu, rlambda, rho;

  Material (int id, std::string style, std::vector<double> param);
  void print() ;

private:
  void set_param();
};
