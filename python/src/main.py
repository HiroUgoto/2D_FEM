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

fp = 3.0
amp = 2.0
printmesh = 2
print("Input frequency(Hz):",fp,"Input amplitude(m/s2):",amp)

## --- EP Set up --- ##
fem.set_ep_initial_state()
# fem.set_rayleigh_damping(fp,10*fp,0.002)

## --- Define input wave --- ##
fsamp = 8000
datasamplerate = 1000  #sampling rate[Hz]
# duration = 5.0/fp + 1.0/fp
# duration = 8.0/fp + 1.0/fp
# duration = 14.0/fp + 1.0/fp
duration = 15

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
h1_h2e_list = np.zeros(fem.output_nelem+1) #h1-h2*e
H1_list = np.zeros(fem.output_nelem+1)
H2_list = np.zeros(fem.output_nelem+1)
L1_list = np.zeros(fem.output_nelem+1)



dispx = open(output_dir+"x.disp","w") 
dispz = open(output_dir+"z.disp","w") 
velx = open(output_dir+"x.vel","w") 
velz = open(output_dir+"z.vel","w") 
accx = open(output_dir+"x.acc","w") 
stressxx = open(output_dir+"stressxx","w")
stresszz = open(output_dir+"stresszz","w")
stressxz = open(output_dir+"stressxz","w")
strainxx = open(output_dir+"xx.str","w")
strainzz = open(output_dir+"zz.str","w")
strainxz = open(output_dir+"xz.str","w")
stressyy = open(output_dir+"stressyy","w")
stresspw = open(output_dir+"stresspw","w")
fL = open(output_dir+"fL.out","w")
psi = open(output_dir+"psi.out","w")
e = open(output_dir+"e.out","w")
n = open(output_dir+"n.out","w")
h = open(output_dir+"h.out","w")
h1_h2e = open(output_dir+"h1-h2e.out","w")
H1 = open(output_dir+"H1.out","w")
H2 = open(output_dir+"H2.out","w")
L1 = open(output_dir+"L1.out","w")



acc0 = np.zeros([2])
vel0 = np.zeros([2])

for it in range(ntim):
    acc0 = np.array([wave_acc[it],0.0])
    vel0 += acc0*dt

    # fem.update_time(acc0,vel0,input_wave=True,self_gravity=True)
    fem.update_time(acc0)
    # fem.update_time(acc0,FD=True)

    if it%(fsamp/datasamplerate) == 0:
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
        writef(e_list, e, tim[it], [element.ep.e for element in fem.output_elements])
        writef(n_list, n, tim[it], [element.ep.n for element in fem.output_elements])
        writef(h_list, h, tim[it], [element.ep.model.h for element in fem.output_elements])
        writef(h1_h2e_list, h1_h2e, tim[it], [(element.ep.model.h1-element.ep.model.h1*element.ep.e) for element in fem.output_elements])
        writef(H1_list, H1, tim[it], [element.ep.model.H1 for element in fem.output_elements])
        writef(H2_list, H2, tim[it], [element.ep.model.H2 for element in fem.output_elements])
        writef(L1_list, L1, tim[it], [element.ep.model.L1 for element in fem.output_elements])
        

        # plot_model.plot_mesh_update(ax,fem,100.)
        print("t=",np.array([it*dt,output_element_stress_xx[printmesh+1],output_element_stress_yy[printmesh+1],output_element_stress_zz[printmesh+1],output_element_stress_zz[printmesh+1]-output_element_pw[printmesh+1],output_element_pw[printmesh+1]]))


dispx.close()
dispz.close()
velx.close()
velz.close()
accx.close()
stressxx.close()
stresszz.close()
stressxz.close()
strainxx.close()
strainzz.close()
strainxz.close()
stressyy.close()
stresspw.close()
fL.close()
psi.close()
e.close()
n.close()
h.close()
h1_h2e.close()
H1.close()
H2.close()
L1.close()


elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")
# plot_model.plot_mesh_update(ax,fem,100.,fin=True)
