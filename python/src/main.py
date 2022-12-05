import matplotlib.pyplot as plt
import numpy as np
import time

import io_data
import input_wave
import plot_model
import datetime
import sys

np.set_printoptions(precision=2,suppress=True)

start = time.time()

## --- Input FEM Mesh --- ##
fem = io_data.input_mesh("input/mesh.in")
outputs = io_data.input_outputs("input/output.in")
output_dir = "result/"

## --- FEM Set up --- ##
fem.set_init()
fem.set_output(outputs)
# plot_model.plot_mesh(fem)

## --- Define input wave --- ##
# Set amplitude #
# vs0,rho0 = 300.0,1700.0  # basement
# vs1,rho1 = 150.0,1700.0  # ground
# h = 10.0    # thickness

# R = (vs1*rho1)/(vs0*rho0)
# omega = 2*np.pi*fp
# c,s = np.cos(omega*h/vs1),np.sin(omega*h/vs1)
# H = 2/np.sqrt(c**2+R**2*s**2)

fp = 2.0
amp = 1.0
print("Input frequency(Hz):",fp,"Input amplitude(m/s2):",amp)

## --- EP Set up --- ##
fem.set_ep_initial_state()
# fem.set_rayleigh_damping(fp,10*fp,0.002)

## --- Define input wave --- ##
fsamp = 10000
# duration = 5.0/fp + 1.0/fp
# duration = 8.0/fp + 1.0/fp
# duration = 14.0/fp + 1.0/fp
duration = 1.0

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
# wave_acc = input_wave.tapered_sin(tim,fp,1.0/fp,duration-1.0/fp,amp)
wave_acc = input_wave.tapered_sin(tim,fp,3.0,duration+3,amp)
# wave_acc = input_wave.simple_sin(tim,fp,amp)
# wave_acc = np.zeros(len(tim))
ntim = len(tim)

# plt.figure()
# plt.plot(tim,wave_acc)
# plt.show()

## --- Prepare time solver --- ##
# ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Iteration ##
def writef(line,outf,tim,data):
    line[0] = tim
    line[1:] = data
    np.savetxt(outf,line[np.newaxis])
    outf.flush()

output_element_stress_xx = np.zeros(fem.output_nelem+1)
output_element_stress_zz = np.zeros(fem.output_nelem+1)
output_element_stress_xz = np.zeros(fem.output_nelem+1)
output_element_stress_yy = np.zeros(fem.output_nelem+1)
output_element_strain_xx = np.zeros(fem.output_nelem+1)
output_element_strain_zz = np.zeros(fem.output_nelem+1)
output_element_strain_xz = np.zeros(fem.output_nelem+1)
output_element_pw = np.zeros(fem.output_nelem+1)

output_accx = np.zeros(fem.output_nnode+1)
output_dispx = np.zeros(fem.output_nnode+1)
output_dispz = np.zeros(fem.output_nnode+1)
output_velx = np.zeros(fem.output_nnode+1)
output_velz = np.zeros(fem.output_nnode+1)

fL_list = np.zeros(fem.output_nelem+1)
psi_list = np.zeros(fem.output_nelem+1)
e_list = np.zeros(fem.output_nelem+1)
n_list = np.zeros(fem.output_nelem+1)
h_list = np.zeros(fem.output_nelem+1)
h12_list = np.zeros(fem.output_nelem+1) #h1-h2*e

with open(output_dir+"disp.x","w") as dispx, \
     open(output_dir+"disp.z","w") as dispz, \
     open(output_dir+"vel.x","w") as velx, \
     open(output_dir+"vel.z","w") as velz, \
     open(output_dir+"acc.x","w") as accx, \
     open(output_dir+"stressxx","w") as stressxx, \
     open(output_dir+"stresszz","w") as stresszz, \
     open(output_dir+"stressxz","w") as stressxz, \
     open(output_dir+"strainxx","w") as strainxx, \
     open(output_dir+"strainzz","w") as strainzz, \
     open(output_dir+"strainxz","w") as strainxz, \
     open(output_dir+"stressyy","w") as stressyy, \
     open(output_dir+"stresspw","w") as stresspw, \
     open(output_dir+"fL","w") as fL, \
     open(output_dir+"psi","w") as psi, \
     open(output_dir+"e","w") as out_e, \
     open(output_dir+"n","w") as out_n, \
     open(output_dir+"h12","w") as h12, \
     open(output_dir+"h","w") as out_h:


    acc0 = np.zeros([2])
    vel0 = np.zeros([2])

    for it in range(ntim):
        acc0 = np.array([wave_acc[it],0.0])
        vel0 += acc0*dt

        # fem.update_time(acc0,vel0,input_wave=True,self_gravity=True)
        fem.update_time(acc0)
        # fem.update_time(acc0,FD=True)

        if it%(fsamp/100) == 0:  #sampling rate 100Hz
            writef(output_accx, accx, tim[it], [node.a[0] for node in fem.output_nodes] + acc0[0])
            writef(output_dispx, dispx, tim[it], [node.u[0] for node in fem.output_nodes])
            writef(output_dispz, dispz, tim[it], [node.u[1] for node in fem.output_nodes])
            writef(output_velx, velx, tim[it], [node.v[0] for node in fem.output_nodes])
            writef(output_velz, velz, tim[it], [node.v[1] for node in fem.output_nodes])
            writef(output_element_stress_xx, stressxx, tim[it], [element.eff_stress[0] for element in fem.output_elements])
            writef(output_element_stress_zz, stresszz, tim[it], [element.eff_stress[1] for element in fem.output_elements])
            writef(output_element_stress_xz, stressxz, tim[it], [element.eff_stress[2] for element in fem.output_elements])
            writef(output_element_strain_xx, strainxx, tim[it], [element.strain[0] for element in fem.output_elements])
            writef(output_element_strain_zz, strainzz, tim[it], [element.strain[1] for element in fem.output_elements])
            writef(output_element_strain_xz, strainxz, tim[it], [element.strain[2] for element in fem.output_elements])
            writef(output_element_stress_yy, stressyy, tim[it], [element.eff_stress_yy for element in fem.output_elements])
            writef(output_element_pw, stresspw, tim[it], [element.excess_pore_pressure for element in fem.output_elements])

            writef(fL_list, fL, tim[it], [element.ep.model.fL for element in fem.output_elements])
            writef(psi_list, psi, tim[it], [element.ep.model.psi for element in fem.output_elements])
            writef(e_list, out_e, tim[it], [element.ep.e for element in fem.output_elements])
            writef(n_list, out_n, tim[it], [element.ep.n for element in fem.output_elements])
            writef(h_list, out_h, tim[it], [element.ep.model.h for element in fem.output_elements])
            writef(h12_list, h12, tim[it], [(element.ep.model.h1-element.ep.model.h1*element.ep.e) for element in fem.output_elements])
            

            # plot_model.plot_mesh_update(ax,fem,100.)
            print("t=",np.array([it*dt,output_element_stress_xx[-1],output_element_stress_yy[-1],output_element_stress_zz[-1],output_element_stress_zz[-1]-output_element_pw[-1],output_element_pw[-1]]))



    elapsed_time = time.time() - start
    print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

    # plot_model.plot_mesh_update(ax,fem,100.,fin=True)


    ## --- Write output file --- ##
    # output_line = np.vstack([element.xnT[:,8] for element in fem.output_elements])
    # np.savetxt(output_dir+"output_element_list.dat",output_line)