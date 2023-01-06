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
  int p = 0;  //標準出力の要素id
  int dsamp = 100;  //出力サンプリングレート

  // ----- FEM Set up ----- //
  fem.set_init();
  fem.set_output(outputs);

  // ----- Define input wave ----- //
  double fp = 3.0;
  double amp = 1.0;

  std::cout << "Input frequency(Hz): " << fp << ", ";
  std::cout << "Input amplitude(m/s2): " << amp << "\n";

  // ----------------------------- //
  fem.set_ep_initial_state();
  fem.set_rayleigh_damping(fp,10*fp,0.02);

  // ----------------------------- //
  size_t fsamp = 20000;
  double duration = 5;

  EV wave_acc;
  auto [tim, dt] = input_wave::linspace(0,duration,(int)(fsamp*duration));
  // wave_acc = input_wave::tapered_sin(tim,fp,1.0/fp,duration-1.0/fp,amp);
  // wave_acc = input_wave::tapered_sin(tim,fp,2.0/fp,duration-1.0/fp,amp);
  wave_acc = input_wave::tapered_sin(tim,fp,1.0,duration+3,amp);
  size_t ntim = tim.size();

  // ----- Prepare time solver ----- //
  fem.update_init(dt);

  // ----- time iteration ----- //
  EV acc0 = EV::Zero(fem.dof);
  EV vel0 = EV::Zero(fem.dof);
  EV dis0 = EV::Zero(fem.dof);

  std::ofstream output_dispx;
  std::ofstream output_dispz;
  std::ofstream output_velx;
  std::ofstream output_velz;
  std::ofstream output_accx;
  std::ofstream output_accz;
  std::ofstream output_strainxx;
  std::ofstream output_strainzz;
  std::ofstream output_strainxz;
  std::ofstream output_stressxx;
  std::ofstream output_stressyy;
  std::ofstream output_stresszz;
  std::ofstream output_stressxz;
  std::ofstream output_pw;
  std::ofstream output_accin;
  std::ofstream output_velin;
  std::ofstream output_dispin;

  // std::ofstream fL;
  // std::ofstream psi;
  // std::ofstream e;
  // std::ofstream n;
  // std::ofstream h;
  // std::ofstream h1_h2e;
  // std::ofstream H1;
  // std::ofstream H2;
  // std::ofstream L1;

  output_dispx.open("result/x.disp");
  output_dispz.open("result/z.disp");
  output_velx.open("result/x.vel");
  output_velz.open("result/z.vel");
  output_accx.open("result/x.acc");
  output_accz.open("result/z.acc");
  output_strainxx.open("result/xx.str");
  output_strainzz.open("result/zz.str");
  output_strainxz.open("result/xz.str");
  output_stressxx.open("result/stressxx");
  output_stressyy.open("result/stressyy");
  output_stresszz.open("result/stresszz");
  output_stressxz.open("result/stressxz");  
  output_pw.open("result/stresspw");
  output_accin.open("result/acc.in");
  output_velin.open("result/vel.in");
  output_dispin.open("result/disp.in");

  // fL.open("result/fL_out");
  // psi.open("result/psi_out");
  // e.open("result/e_out");
  // n.open("result/n_out");
  // h.open("result/h_out");
  // h1_h2e.open("result/h1_h2e_out");
  // H1.open("result/H1_out");
  // H2.open("result/H2_out");
  // L1.open("result/L1_out");

  for (size_t it = 0 ; it < ntim ; it++) {
    acc0[0] = wave_acc[it];
    vel0[0] += wave_acc[it]*dt;


    fem.update_time_MD(acc0);
    // fem.update_time_FD(acc0);
    // fem.update_time_input_MD(vel0);
    // fem.update_time_input_FD(vel0);

  // --- Write output file --- //
    if (it%(fsamp/dsamp) == 0){
      output_dispx << tim(it);
      output_dispz << tim(it);
      output_velx << tim(it);
      output_velz << tim(it);
      output_accx << tim(it);
      output_accz << tim(it);
      output_strainxx << tim(it);
      output_strainzz << tim(it);
      output_strainxz << tim(it);
      output_stressxx << tim(it);
      output_stresszz << tim(it);
      output_stressxz << tim(it);
      output_stressyy << tim(it);
      output_pw << tim(it);
      output_accin << acc0[0] << std::endl;
      output_velin << vel0[0] << std::endl;
      output_dispin << dis0[0] << std::endl;

      for (size_t i = 0 ; i < fem.output_nnode ; i++) {
        Node* node_p = fem.output_nodes_p[i];
        output_dispx << " " << node_p->u(0);
        output_dispz << " " << node_p->u(1);
        output_velx << " " << node_p->v(0);
        output_velz << " " << node_p->v(1);
        output_accx << " " << node_p->a(0) + acc0[0];   //慣性力入力
        // output_accx << " " << node_p->a(0);   //input要素
        output_accz << " " << node_p->a(1);
      }

      for (size_t i = 0 ; i < fem.output_nelem ; i++) {
        Element* element_p = fem.output_elements_p[i];
        output_strainxx << " " << element_p->strain(0);
        output_strainzz << " " << element_p->strain(1);
        output_strainxz << " " << element_p->strain(2);
        output_stressxx << " " << element_p->eff_stress(0);
        output_stresszz << " " << element_p->eff_stress(1);
        output_stressyy << " " << element_p->eff_stress_yy;
        output_stressxz << " " << element_p->eff_stress(2);
        output_pw << " " << element_p->excess_pore_pressure;
      }
      
      output_dispx << std::endl;
      output_dispz << std::endl;
      output_velx << std::endl;
      output_velz << std::endl;
      output_accx << std::endl;
      output_accz << std::endl;
      output_strainxx << std::endl;
      output_strainzz << std::endl;
      output_strainxz << std::endl;
      output_stressxx << std::endl;
      output_stresszz << std::endl;
      output_stressyy << std::endl;
      output_stressxz << std::endl;
      output_pw << std::endl;

      if (it%(fsamp/dsamp*1) == 0){
        Element* out_element_p = &fem.elements[p];
        // std::cout << " t=" << it*dt << " e_id:" << p << " " << out_element_p->stress(0) << " " << out_element_p->stress(1) << std::endl;
        std::cout << " t=" << it*dt << " e_id:" << p << " " << out_element_p->eff_stress(0) << " " << out_element_p->eff_stress_yy << " " << out_element_p->eff_stress(1) << " " << out_element_p->excess_pore_pressure << std::endl;
      }
    }
  }

  output_dispx.close();
  output_dispz.close();
  output_velx.close();
  output_velz.close();
  output_accx.close();
  output_accz.close();

  output_strainxx.close();
  output_strainzz.close();
  output_strainxz.close();

  output_stressxx.close();
  output_stressyy.close();
  output_stresszz.close();
  output_stressxz.close();
  output_pw.close();

  output_accin.close();
  output_velin.close();
  output_dispin.close();

  clock_t end = clock();
  std::cout << "elapsed_time: " << (double)(end - start) / CLOCKS_PER_SEC << "[sec]\n";
}
