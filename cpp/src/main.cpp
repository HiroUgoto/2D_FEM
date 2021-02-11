#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "element.h"
#include "fem.h"
#include "io_data.h"
#include "input_wave.h"

int main() {

  clock_t start = clock();

  // ----- Input FEM Mesh ----- //
  Fem fem = io_data::input_mesh("input/mesh.in");
  auto outputs = io_data::input_outputs("input/output.in");

  // ----- FEM Set up ----- //
  fem.set_init();
  fem.set_output(outputs);

  // ----- Define input wave ----- //
  size_t fsamp = 5000;
  double fp = 0.2;
  double duration = 1.0/fp;

  Eigen::VectorXd tim, wave_acc;
  double dt;
  std::tie(tim, dt) = input_wave::linspace(0,duration,(int)(fsamp*duration));
  wave_acc = input_wave::simple_sin(tim,fp,0.1);
  size_t ntim = tim.size();

  // ----- Static deformation ----- //


  clock_t end = clock();
  std::cout << "elapsed_time: " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n";

}
