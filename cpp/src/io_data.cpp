#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "fem.h"
#include "io_data.h"

// ------------------------------------------------------------------- //
Fem io_data::input_mesh (const std::string mesh_file) {
    size_t nnode, nelem, nmaterial, dof;
    std::string line;

    // Open file //
    std::ifstream f(mesh_file);

    // Read header //
    std::getline(f, line);
    std::istringstream iss(line);
    iss >> nnode >> nelem >> nmaterial >> dof;

    // Read nodes //
    std::vector<Node> nodes;
    for (size_t inode = 0 ; inode < nnode ; inode++) {
      size_t id;
      std::vector<double> xyz(2);
      std::vector<size_t> freedom(dof);

      std::getline(f, line);
      // std::cout << line + "\n";

      std::istringstream iss(line);
      iss >> id;
      for (size_t i = 0 ; i < 2 ; i++) {
        iss >> xyz.at(i);
      }
      for (size_t i = 0 ; i < dof ; i++) {
        iss >> freedom.at(i);
      }

      Node node(id,xyz,freedom);
      nodes.push_back(node);
    }

    // Read elements //
    std::vector<Element> elements;
    for (size_t ielem = 0 ; ielem < nelem ; ielem++) {
      size_t id;
      int material_id;
      std::string style;
      std::vector<size_t> inode;

      std::getline(f, line);
      // std::cout << line + "\n";

      std::istringstream iss(line);
      iss >> id >> style >> material_id ;
      while(!iss.eof()) {
        size_t in;
        iss >> in;
        inode.push_back(in);
      }

      Element element(id,style,material_id,inode);
      elements.push_back(element);
    }

    // Read materials //
    std::vector<Material> materials;
    for (size_t imaterial = 0 ; imaterial < nmaterial ; imaterial++) {
      size_t id;
      std::string style;
      std::vector<double> param;

      std::getline(f, line);
      // std::cout << line + "\n";

      std::istringstream iss(line);
      iss >> id >> style ;
      while(!iss.eof()) {
        double ip;
        iss >> ip;
        param.push_back(ip);
      }

      Material material(id,style,param);
      materials.push_back(material);
    }

    Fem fem(dof,nodes,elements,materials);
    return fem;
  }

// ------------------------------------------------------------------- //
std::tuple<std::vector<size_t>, std::vector<size_t>>
  io_data::input_outputs (const std::string output_file) {
    size_t nnode, nelem;
    std::string line;

    // Open file //
    std::ifstream f(output_file);

    // Read header //
    std::getline(f, line);
    std::istringstream iss(line);
    iss >> nnode >> nelem;

    // Read nodes //
    std::vector<size_t> nodes;
    for (size_t inode = 0 ; inode < nnode ; inode++) {
      size_t id;

      std::getline(f, line);
      // std::cout << line + "\n";
      std::istringstream iss(line);
      iss >> id;

      nodes.push_back(id);
    }

    std::vector<size_t> elements;
    for (size_t ielem = 0 ; ielem < nelem ; ielem++) {
      size_t id;

      std::getline(f, line);
      // std::cout << line + "\n";
      std::istringstream iss(line);
      iss >> id;

      elements.push_back(id);
    }

    return {nodes, elements};
  }
