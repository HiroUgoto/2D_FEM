#include "all.h"
#include <Eigen/Core>
#include "node.h"
#include "material.h"
#include "element_style.h"
#include "elasto_plastic.h"
#include "ep_model.h"
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
  double vs0 = 300;
  double rho0 = 1700;
  double vs1 = 150;
  double rho1 = 1700;
  double h = 10;

  double fp = vs1/(4*h);
  double R = (vs1*rho1)/(vs0*rho0);
  double omega = 2*M_PI*fp;
  double H = 2/sqrt(pow(cos(omega*h/vs1),2) + R*R*pow(sin(omega*h/vs1),2));
  double amp = 0.3*9.8 / H;

  std::cout << "Input frequency(Hz): " << fp << ", ";
  std::cout << "Input amplitude(m/s2): " << amp << "\n";

  // ----------------------------- //
  fem.set_ep_initial_state();
  fem.set_rayleigh_damping(fp,10*fp,0.002);

  // ----------------------------- //
  size_t fsamp = 15000;
  amp = 0.75 ;
  fp = 3.75 ;
  // double duration = 5.0/fp + 1.0/fp;
  // double duration = 8.0/fp + 1.0/fp;
  double duration = 14.0/fp + 1.0/fp;

  EV wave_acc;
  auto [tim, dt] = input_wave::linspace(0,duration,(int)(fsamp*duration));
  // wave_acc = input_wave::tapered_sin(tim,fp,1.0/fp,duration-1.0/fp,amp);
  // wave_acc = input_wave::tapered_sin(tim,fp,2.0/fp,duration-1.0/fp,amp);
  wave_acc = input_wave::tapered_sin(tim,fp,3.0/fp,duration-1.0/fp,amp);
  size_t ntim = tim.size();

  std::ofstream f0(output_dir + "input.acc");
  for (size_t it = 0 ; it < ntim ; it++) {
    f0 << tim(it) ;
    f0 << " " << wave_acc(it) ;
    f0 << "\n";
  }
  f0.close();

  // ----- Prepare time solver ----- //
  fem.update_init(dt);

  EM output_dispx = EM::Zero(ntim,fem.output_nnode);
  EM output_dispz = EM::Zero(ntim,fem.output_nnode);
  EM output_accx = EM::Zero(ntim,fem.output_nnode);
  EM output_accz = EM::Zero(ntim,fem.output_nnode);

  EM output_element_stress_xx = EM::Zero(ntim,fem.output_nelem);
  EM output_element_stress_zz = EM::Zero(ntim,fem.output_nelem);
  EM output_element_stress_xz = EM::Zero(ntim,fem.output_nelem);
  EM output_element_stress_yy = EM::Zero(ntim,fem.output_nelem);

  // ----- time iteration ----- //
  EV acc0 = EV::Zero(fem.dof);
  EV vel0 = EV::Zero(fem.dof);

  // for (size_t it = 0 ; it < 1000 ; it++) {
  for (size_t it = 0 ; it < ntim ; it++) {
    acc0[0] = wave_acc[it];
    vel0[0] += wave_acc[it]*dt;

    // fem.update_time_input_MD(vel0);
    fem.update_time_input_MD_gravity(vel0);
    // fem.update_time_MD_gravity(acc0);

    for (size_t i = 0 ; i < fem.output_nnode ; i++) {
      Node* node_p = fem.output_nodes_p[i];
      output_dispx(it,i) = node_p->u(0);
      output_dispz(it,i) = node_p->u(1);
      output_accx(it,i) = node_p->a(0) + acc0[0];
      output_accz(it,i) = node_p->a(1);
    }

    for (size_t i = 0 ; i < fem.output_nelem ; i++) {
      Element* element_p = fem.output_elements_p[i];
      output_element_stress_xx(it,i) = element_p->stress(0);
      output_element_stress_zz(it,i) = element_p->stress(1);
      output_element_stress_xz(it,i) = element_p->stress(2);
      output_element_stress_yy(it,i) = element_p->stress_yy;
    }

    // std::cout << fem.elements[0].stress(0) << " ";
    // std::cout << fem.elements[0].stress(1) << " ";
    // std::cout << fem.elements[0].stress(2) << " ";
    // std::cout << fem.elements[0].stress_yy << std::endl;

    if (it%100 == 0) {
      std::cout << it << " t= " << it*dt << " ";
      std::cout << output_accx(it,0) << " ";
      std::cout << output_element_stress_xx(it,1) << " ";
      std::cout << output_element_stress_zz(it,1) << " ";
      std::cout << output_element_stress_yy(it,1) << std::endl;
    }
  }

  clock_t end = clock();
  std::cout << "elapsed_time: " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n";

  // --- Write output file --- //
  std::ofstream fa(output_dir + "result.acc");
  std::ofstream fd(output_dir + "result.disp");
  for (size_t it = 0 ; it < ntim ; it++) {
    fa << tim(it) ;
    fa << " " << output_accx(it,0);
    fa << "\n";

    fd << tim(it) ;
    fd << " " << output_dispx(it,0);
    fd << " " << output_dispx(it,int(fem.output_nnode/2));
    fd << "\n";
  }
  fa.close();
  fd.close();

  std::ofstream f(output_dir + "output_element_list.dat");
  for (size_t i = 0 ; i < fem.output_nelem ; i++) {
    Element* element_p = fem.output_elements_p[i];
    f << element_p->xnT(0,8) << " " << element_p->xnT(1,8) << "\n";
  }
  f.close();

  std::ofstream fxx(output_dir + "output_element.stress_xx");
  std::ofstream fzz(output_dir + "output_element.stress_zz");
  std::ofstream fxz(output_dir + "output_element.stress_xz");
  std::ofstream fyy(output_dir + "output_element.stress_yy");
  for (size_t it = 0 ; it < ntim ; it++) {
    fxx << tim(it) ;
    fzz << tim(it) ;
    fxz << tim(it) ;
    fyy << tim(it) ;
    for (size_t i = 0 ; i < fem.output_nelem ; i++) {
      fxx << " " << output_element_stress_xx(it,i);
      fzz << " " << output_element_stress_zz(it,i);
      fxz << " " << output_element_stress_xz(it,i);
      fyy << " " << output_element_stress_yy(it,i);
    }
    fxx << "\n";
    fzz << "\n";
    fxz << "\n";
    fyy << "\n";
  }
  fxx.close();
  fzz.close();
  fxz.close();
  fyy.close();
}
