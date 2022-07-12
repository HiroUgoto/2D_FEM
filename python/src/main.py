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
fsamp = 1000
fp = 2.5
duration = 4.0/fp

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
slip = input_wave.smoothed_ramp(tim,fp=fp,tp=1.0/fp,amp=1.0)
# slip_rate = input_wave.ricker(tim,fp=fp,tp=1.0/fp,amp=1.0)
ntim = len(tim)


# plt.figure()
# plt.plot(tim,slip)
# plt.show()

## --- Prepare time solver --- ##
ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Iteration ##
output_dispx = np.zeros((ntim,fem.output_nnode))
output_dispz = np.zeros((ntim,fem.output_nnode))

output_velx = np.zeros((ntim,fem.output_nnode))
output_velz = np.zeros((ntim,fem.output_nnode))

output_accx = np.zeros((ntim,fem.output_nnode))
output_accz = np.zeros((ntim,fem.output_nnode))

# slip0 = 0.0
for it in range(len(tim)):
    slip0 = slip[it]
    # slip_rate0 = slip_rate[it]
    # slip0 += slip_rate0*dt

    fem.update_time_slip(slip0)

    output_dispx[it,:] = [node.u[0] for node in fem.output_nodes]
    output_dispz[it,:] = [node.u[1] for node in fem.output_nodes]

    output_velx[it,:] = [node.v[0] for node in fem.output_nodes]
    output_velz[it,:] = [node.v[1] for node in fem.output_nodes]

    output_accx[it,:] = [node.a[0] for node in fem.output_nodes]
    output_accz[it,:] = [node.a[1] for node in fem.output_nodes]

    if it%20 == 0:
        plot_model.plot_mesh_update(ax,fem,1.)
        print(it,"t=",it*dt,output_dispx[it,int(fem.output_nnode//2)])

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# plot_model.plot_mesh_update(ax,fem,10.,fin=True)

## --- Write output file --- ##
output_line = np.vstack([tim,output_dispx[:,0],output_dispx[:,int(fem.output_nnode//2)]]).T
np.savetxt(output_dir+"output_x.disp",output_line)

output_line = np.vstack([tim,output_velx[:,0],output_velx[:,int(fem.output_nnode//2)]]).T
np.savetxt(output_dir+"output_x.vel",output_line)

output_line = np.vstack([tim,output_accx[:,0],output_accx[:,int(fem.output_nnode//2)]]).T
np.savetxt(output_dir+"output_x.acc",output_line)

## Output result ##
plt.figure()
plt.plot(tim,wave_acc,c='k')
plt.plot(tim,output_accx[:,0],c='r')
plt.show()
