import matplotlib.pyplot as plt
import numpy as np
import time

import io_data
import input_wave
import plot_model

start = time.time()

## --- Input FEM Mesh --- ##
fem = io_data.input_mesh("input/mesh.in")
outputs = io_data.input_outputs("input/output.in")

## --- FEM Set up --- ##
fem.set_init()
fem.set_output(outputs)
# plot_model.plot_mesh(fem)


## --- Define input wave --- ##
fsamp = 2500
duration = 5.0
theta = 10 * np.pi/180.0

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
horizontal_acc = input_wave.smoothed_ramp(tim,fp=5.0,tp=1/5.,amp=9.8*np.sin(theta))
vertical_acc = -input_wave.smoothed_ramp(tim,fp=5.0,tp=1/5.,amp=9.8*(1-np.cos(theta)))

# wave_acc = input_wave.tapered_sin(tim,fp=5.0,taper=0.2,duration=1.0,amp=1.0)
# wave_vel = np.cumsum(wave_acc) * dt
ntim = len(tim)

ax = plot_model.plot_mesh_update_init()
## --- Static deformation --- ##
fem.self_gravity()
plot_model.plot_mesh_update(ax,fem,1.)

## --- Prepare time solver --- ##
fem.update_init(dt)

## Iteration ##
output_disp = np.zeros((ntim,fem.output_nnode))
output_vel = np.zeros((ntim,fem.output_nnode))
output_strain = np.zeros((ntim,fem.output_nelem))

for it in range(len(tim)):
    acc0 = np.array([horizontal_acc[it],vertical_acc[it]])

    fem.update_time(acc0,FD=True)
    # fem.update_time(acc0)

    output_disp[it,:] = [node.u[1]-node.u0[1] for node in fem.output_nodes]
    output_vel[it,:] = [node.v[0] for node in fem.output_nodes]
    output_strain[it,:] = [element.strain[0] for element in fem.output_elements]

    if it%10 == 0:
        plot_model.plot_mesh_update(ax,fem,1.)
        print(it,"t=",it*dt,output_disp[it,:]," ans =",10*np.tan(theta))


plot_model.plot_mesh_update(ax,fem,1.,fin=True)

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# ## Output result ##
# plt.figure()
# plt.plot(tim,wave_vel)
# plt.plot(tim,output_vel[:,0])
# plt.show()
