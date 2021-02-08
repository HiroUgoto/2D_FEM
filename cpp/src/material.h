class Material {
public:
  int id;
  std::string style;
  std::vector<double> param;

  Material (int id, std::string style, std::vector<double> param);
  void print() ;
};
