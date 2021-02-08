class Element {
public:
  int id;
  int material_id;
  std::vector<int> inode;
  std::string style;
  double gravity;

  Element (int id, std::string style, int material_id, std::vector<int> inode);
  void print() ;
};
