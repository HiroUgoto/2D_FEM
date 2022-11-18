import matplotlib.pyplot as plt
import numpy as np
import time

import io_data
import input_wave
import plot_model
import datetime
import sys,os

from ep_model.DL1d import constitution as DL1d

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
# fp = 0.5  # Hz
fp = 1  # Hz
amp = 1.0 # m/s
print("Input frequency(Hz):",fp,"Input amplitude(m/s2):",amp)

## --- EP Set up --- ##
fem.set_ep_initial_state()
# fem.set_rayleigh_damping(fp,10*fp,0.002)

## --- Define input wave --- ##
fsamp = 2000
duration = 10.0/fp + 1.0/fp

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
wave_acc = input_wave.tapered_sin(tim,fp,3.0/fp,duration-1.0/fp,amp)
ntim = len(tim)

# plt.figure()
# plt.plot(tim,wave_acc)
# plt.show()

## --- Prepare time solver --- ##
ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Iteration ##
output_accx = np.zeros((ntim,fem.output_nnode))
output_element_strain_xz = np.zeros((ntim,fem.output_nelem))
output_element_stress_xz = np.zeros((ntim,fem.output_nelem))

acc0 = np.array([0.0,0.0])
vel0 = np.array([0.0,0.0])

for it in range(ntim):
    acc0 = np.array([wave_acc[it],0.0])
    vel0 += acc0*dt

    fem.update_time(acc0,vel0,input_wave=True,self_gravity=True)

    output_accx[it,:] = [node.a[0] for node in fem.output_nodes]

    output_element_strain_xz[it,:] = [element.strain[2] for element in fem.output_elements]
    output_element_stress_xz[it,:] = [element.stress[2] for element in fem.output_elements]

    if it%50 == 0:
        plot_model.plot_mesh_update(ax,fem,10.,margin=2)
        print(it,"t=",tim[it],output_accx[it,0],output_element_strain_xz[it,0],output_element_stress_xz[it,0])

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# plot_model.plot_mesh_update(ax,fem,10.,fin=True,margin=0.5)

## --- Write output file --- ##
output_line = np.vstack([tim[:ntim],wave_acc[:ntim]]).T
np.savetxt(output_dir+"input.acc",output_line)

output_line = np.vstack([tim[:ntim],output_accx[:,0]]).T
np.savetxt(output_dir+"result.acc",output_line)

for ielem in range(fem.output_nelem):
    output_line = np.vstack([tim[:ntim],output_element_strain_xz[:,ielem]])
np.savetxt(output_dir+"output_element.strain_xz",output_line.T)

for ielem in range(fem.output_nelem):
    output_line = np.vstack([tim[:ntim],output_element_stress_xz[:,ielem]])
np.savetxt(output_dir+"output_element.stress_xz",output_line.T)

for i,element in enumerate(fem.output_elements):
    fname = "result/fig/gamma_tau_{:0>2}.png".format(i)
    element.ep.model.plot(fname)
