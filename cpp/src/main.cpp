#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "fem.h"
#include "io_data.h"
#include "input_wave.h"

using EV = Eigen::VectorXd ;
using EM = Eigen::MatrixXd ;

int main() {

  clock_t start = clock();

  // ----- Input FEM Mesh ----- //
  Fem fem = io_data::input_mesh("input/mesh.in");
  auto outputs = io_data::input_outputs("input/output.in");
  std::string output_dir = "result/";

  // ----- FEM Set up ----- //
  fem.set_init();
  fem.set_output(outputs);

  // ----- Define input wave ----- //
  size_t fsamp = 10000;
  double fp = 0.196;
  double duration = 8.0/fp;

  EV wave_acc;
  auto [tim, dt] = input_wave::linspace(0,duration,(int)(fsamp*duration));
  wave_acc = input_wave::simple_sin(tim,fp,0.01);
  size_t ntim = tim.size();

  // std::ofstream f0(output_dir + "input.acc");
  // for (size_t it = 0 ; it < ntim ; it++) {
  //   f0 << tim(it) ;
  //   f0 << " " << wave_acc(it) ;
  //   f0 << "\n";
  // }
  // f0.close();
  // exit(1);

  // ----- Static deformation ----- //
  fem.self_gravity();

  // ----- Prepare time solver ----- //
  fem.update_init(dt);

  EM output_dispx = EM::Zero(ntim,fem.output_nnode);
  EM output_dispz = EM::Zero(ntim,fem.output_nnode);
  EM output_velx = EM::Zero(ntim,fem.output_nnode);
  EM output_velz = EM::Zero(ntim,fem.output_nnode);

  // ----- time iteration ----- //
  EV acc0 = EV::Zero(fem.dof);
  EV vel0 = EV::Zero(fem.dof);

  for (size_t it = 0 ; it < ntim ; it++) {
    acc0[0] = wave_acc[it];

    // fem.update_time_MD(acc0);
    fem.update_time_FD(acc0);

    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      Node* node_p = fem.output_nodes_p[i];
      output_dispx(it,i) = node_p->u(0) - node_p->u0(0);
      output_dispz(it,i) = node_p->u(1) - node_p->u0(1);
      output_velx(it,i) = node_p->v(0);
      output_velz(it,i) = node_p->v(1);
    }

    if (it%500 == 0) {
      std::cout << it << " t= " << it*dt << " ";
      std::cout << output_dispz(it,0) << "\n";
    }
  }

  clock_t end = clock();
  std::cout << "elapsed_time: " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n";

  // --- Write output file --- //
  std::ofstream f(output_dir + "z0_vs10_test.disp");
  for (size_t it = 0 ; it < ntim ; it++) {
    f << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      f << " " << output_dispz(it,i);
    }
    f << "\n";
  }
  f.close();

}
