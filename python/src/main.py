import matplotlib.pyplot as plt
import numpy as np
import time

import io_data
import input_wave
import plot_model
import datetime
import sys

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
# # Set amplitude #
# vs0,rho0 = 300.0,1700.0  # basement
# vs1,rho1 = 150.0,1700.0  # ground
# h = 10.0    # thickness
# fp = vs1/(4*h)
#
# R = (vs1*rho1)/(vs0*rho0)
# omega = 2*np.pi*fp
# c,s = np.cos(omega*h/vs1),np.sin(omega*h/vs1)
# H = 2/np.sqrt(c**2+R**2*s**2)
#
# amp = 0.3*9.8 / H
fp = 2.0
amp = 1.0
print("Input frequency(Hz):",fp,", Input amplitude(m/s2):",amp)

## --- EP Set up --- ##
fem.set_ep_initial_state()
# fem.set_rayleigh_damping(fp,10*fp,0.002)

## --- Define input wave --- ##
fsamp = 10000
duration = 5.0/fp + 1.0/fp

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
wave_acc = input_wave.tapered_sin(tim,fp,1.0/fp,duration-1.0/fp,amp)
ntim = len(tim)

# plt.figure()
# plt.plot(tim,wave_acc)
# plt.show()

## --- Prepare time solver --- ##
ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Iteration ##
output_element_stress_xx = np.zeros((ntim,fem.output_nelem))
output_element_stress_zz = np.zeros((ntim,fem.output_nelem))
output_element_stress_xz = np.zeros((ntim,fem.output_nelem))
output_element_stress_yy = np.zeros((ntim,fem.output_nelem))
output_element_pore_pressure = np.zeros((ntim,fem.output_nelem))

output_element_strain_xz = np.zeros((ntim,fem.output_nelem))

output_accx = np.zeros((ntim,fem.output_nnode))
output_dispx = np.zeros((ntim,fem.output_nnode))
output_dispz = np.zeros((ntim,fem.output_nnode))

acc0 = np.array([0.0,0.0])
vel0 = np.array([0.0,0.0])


for it in range(ntim):
    acc0 = np.array([wave_acc[it],0.0])
    vel0 += acc0*dt

    fem.update_time(acc0,vel0,input_wave=True)

    output_accx[it,:] = [node.a[0] for node in fem.output_nodes] + acc0[0]
    output_dispx[it,:] = [node.u[0] for node in fem.output_nodes]
    output_dispz[it,:] = [node.u[1] for node in fem.output_nodes]

    output_element_stress_xx[it,:] = [element.stress[0] for element in fem.output_elements]
    output_element_stress_zz[it,:] = [element.stress[1] for element in fem.output_elements]
    output_element_stress_xz[it,:] = [element.stress[2] for element in fem.output_elements]
    output_element_stress_yy[it,:] = [element.stress_yy for element in fem.output_elements]
    output_element_pore_pressure[it,:] = [element.pore_pressure for element in fem.output_elements]

    output_element_strain_xz[it,:] = [element.strain[2] for element in fem.output_elements]

    if it%50 == 0:
        plot_model.plot_mesh_update(ax,fem,50.)
        # print(it,"t=",it*dt,output_accx[it,0],output_element_stress_xx[it,0],output_element_stress_zz[it,0],output_element_stress_yy[it,0])
        # print(it,"t=",it*dt,output_accx[it,0],output_element_strain_xz[it,0],output_element_stress_xz[it,0])
        print(it,"t=",it*dt,output_dispx[it,0],output_element_stress_xx[it,0],output_element_stress_zz[it,0],output_element_pore_pressure[it,0])


elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# plot_model.plot_mesh_update(ax,fem,100.,fin=True)

## --- Write output file --- ##
output_line = np.vstack([tim[:ntim],output_dispx[:,0],output_dispx[:,int(fem.output_nnode//2)]]).T
np.savetxt(output_dir+"result.disp",output_line)

output_line = np.vstack([tim[:ntim],output_accx[:,0]]).T
np.savetxt(output_dir+"result.acc",output_line)

#
output_line = np.vstack([element.xnT[:,8] for element in fem.output_elements])
np.savetxt(output_dir+"output_element_list.dat",output_line)

output_line = tim[:ntim]
for ielem in range(fem.output_nelem):
    output_line = np.vstack([output_line,output_element_stress_xx[:,ielem]])
np.savetxt(output_dir+"output_element.stress_xx",output_line.T)

output_line = tim[:ntim]
for ielem in range(fem.output_nelem):
    output_line = np.vstack([output_line,output_element_stress_zz[:,ielem]])
np.savetxt(output_dir+"output_element.stress_zz",output_line.T)

output_line = tim[:ntim]
for ielem in range(fem.output_nelem):
    output_line = np.vstack([output_line,output_element_stress_xz[:,ielem]])
np.savetxt(output_dir+"output_element.stress_xz",output_line.T)

output_line = tim[:ntim]
for ielem in range(fem.output_nelem):
    output_line = np.vstack([output_line,output_element_pore_pressure[:,ielem]])
np.savetxt(output_dir+"output_element.pore_pressure",output_line.T)

output_line = tim[:ntim]
for ielem in range(fem.output_nelem):
    output_line = np.vstack([output_line,output_element_strain_xz[:,ielem]])
np.savetxt(output_dir+"output_element.strain_xz",output_line.T)
