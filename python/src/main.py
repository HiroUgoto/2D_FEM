import matplotlib.pyplot as plt
import numpy as np
import time

import io_data
import input_wave
import plot_model
import datetime


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
fsamp = 10000
fp = 1.0/0.264
duration = 14.0/fp + 1.0/fp

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
wave_acc = input_wave.tapered_sin(tim,fp,3.0/fp,14.0/fp,3.0)
ntim = len(tim)

# plt.figure()
# plt.plot(tim,wave_acc)
# plt.show()

## --- Prepare time solver --- ##
ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Iteration ##
num_spring = len(fem.spring_elements)
output_spring_forcev = np.zeros((ntim,num_spring))
output_spring_forceh = np.zeros((ntim,num_spring))

output_dispx = np.zeros((ntim,fem.output_nnode))
output_dispz = np.zeros((ntim,fem.output_nnode))

acc0 = np.array([0.0,0.0])
vel0 = np.array([0.0,0.0])

for it in range(ntim):
    acc0 = np.array([wave_acc[it],0.0])
    vel0 += acc0*dt

    fem.update_time(acc0,vel0,input_wave=True)

    output_dispx[it,:] = [node.u[0] for node in fem.output_nodes]
    output_dispz[it,:] = [node.u[1] for node in fem.output_nodes]

    output_spring_forcev[it,:] = [element.f[0] for element in fem.spring_elements]
    output_spring_forceh[it,:] = [element.f[1] for element in fem.spring_elements]

    if it%500 == 0:
        plot_model.plot_mesh_update(ax,fem,50.)
        print(it,"t=",it*dt,output_dispx[it,5])


elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# plot_model.plot_mesh_update(ax,fem,50.,fin=True)

## --- Write output file --- ##
output_line = np.vstack([tim[:ntim],output_dispx[:,5]]).T
np.savetxt(output_dir+"result_spring.disp",output_line)
#
#
output_line = np.vstack([element.xnT[:,0] for element in fem.spring_elements])
np.savetxt(output_dir+"spring_list.dat",output_line)

output_line = np.vstack([tim[:ntim],output_spring_forcev[:,:].T]).T
np.savetxt(output_dir+"result_spring.forcev",output_line)

## Output result ##
# plt.figure()
# plt.plot(tim,output_dispx[:,5])
# plt.show()
