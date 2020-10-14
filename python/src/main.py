import matplotlib.pyplot as plt
import numpy as np

import io_data
import input_wave
import plot_model

## --- Input FEM Mesh --- ##
fem = io_data.input_mesh("input/mesh.in")
outputs = io_data.input_outputs("input/output.in")

## --- FEM Set up --- ##
fem.set_init()
fem.set_output(outputs)
# plot_model.plot_mesh(fem)

## --- Define input wave --- ##
fsamp = 5000
duration = 0.50

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
wave_acc = input_wave.ricker(tim,0.20,5.0,1.0)
wave_vel = np.cumsum(wave_acc) * dt
ntim = len(tim)

## --- Prepare time solver --- ##
fem.update_init(dt)
ax = plot_model.plot_mesh_update_init()

## Iteration ##
output_vel = np.zeros((ntim,fem.output_nnode))
output_strain = np.zeros((ntim,fem.output_nelem))

for it in range(len(tim)):
    acc0 = np.array([0.0,0.0])
    vel0 = np.array([wave_vel[it],0.0])

    fem.update_time(acc0,vel0,input_wave=True)

    output_vel[it,:] = [node.v[0] for node in fem.output_nodes]
    output_strain[it,:] = [element.strain[0] for element in fem.output_elements]

    if it%50 == 0:
        plot_model.plot_mesh_update(ax,fem,500.)
        print(it,output_vel[it,0],output_strain[it,0])

## Output result ##
plt.figure()
plt.plot(tim,wave_vel)
plt.plot(tim,output_vel[:,0])
plt.show()
