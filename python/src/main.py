import matplotlib.pyplot as plt
import numpy as np
import time

import io_data
import input_wave
import source
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
# exit()

## --- Define EQ source --- ##
fsamp = 100

fp = 1
duration = 12

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
# slip_rate = input_wave.ricker(tim,fp,tp=1.0/fp,amp=1.0)
slip = input_wave.smoothed_ramp(tim,fp,tp=1.0/fp,amp=1.0)
ntim = len(tim)

dip = 30.0  # degree
width = 100.0   # m
sx = 1500   # m
sz = 3000   # m
sources = source.set_source(fem.dof,fem.elements,dip,width,sx,sz)

# plt.figure()
# plt.plot(tim,slip)
# plt.show()
# exit()

## --- Prepare time solver --- ##
ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Iteration ##
output_disp = np.zeros((ntim,fem.output_nnode))
output_vel = np.zeros((ntim,fem.output_nnode))
output_acc = np.zeros((ntim,fem.output_nnode))

slip0 = 0.0
for it in range(len(tim)):
    # slip0 += slip_rate[it]*dt
    slip0 = slip[it]

    fem.update_time_source(sources,slip0)

    output_disp[it,:] = [node.u[0] for node in fem.output_nodes]
    output_vel[it,:] = [node.v[0] for node in fem.output_nodes]
    output_acc[it,:] = [node.a[0] for node in fem.output_nodes]

    if it%10 == 0:
        plot_model.plot_mesh_update(ax,fem,200.)
        print(it,"t=",it*dt,output_disp[it,0])

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# plot_model.plot_mesh_update(ax,fem,10.,fin=True)

## --- Write output file --- ##
output_line = np.vstack([tim,output_disp.T]).T
np.savetxt(output_dir+"output.disp",output_line)

output_line = np.vstack([tim,output_vel.T]).T
np.savetxt(output_dir+"output.vel",output_line)

output_line = np.vstack([tim,output_acc.T]).T
np.savetxt(output_dir+"output.acc",output_line)

## Output result ##
plt.figure()
# plt.plot(tim,slip_rate,c='k')
plt.plot(tim,output_vel[:,0],c='k')
plt.plot(tim,output_vel[:,15],c='r')
plt.plot(tim,output_vel[:,30],c='b')
plt.show()
