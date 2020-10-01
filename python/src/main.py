import matplotlib.pyplot as plt
import numpy as np
import time

import io_data
import input_wave
import plot_model

start = time.time()

## --- Input FEM Mesh --- ##
fem = io_data.input_mesh("input/mesh.in")

## --- FEM Set up --- ##
fem.set_init()
# plot_model.plot_mesh(fem)

## --- Define input wave --- ##
fsamp = 5000
duration = 0.50

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
wave_acc = input_wave.ricker(tim,0.10,10.0,1.0)
wave_vel = np.cumsum(wave_acc) * dt

## --- Prepare time solver --- ##
fem.update_init(dt)
ax = plot_model.plot_mesh_update_init()

## Iteration ##
output_vel = np.zeros_like(wave_vel)
for it in range(len(tim)):
    acc0 = np.array([0.0,0.0])
    vel0 = np.array([wave_vel[it],0.0])

    fem.update_time(acc0,vel0,input_wave=True)
    output_vel[it] = np.copy(fem.nodes[20].v[0])

    if it%50 == 0:
        plot_model.plot_mesh_update(ax,fem,5000.)
        print(it,fem.nodes[20].u)

## Check process time ##
process_time = time.time() - start
print("process time: ", process_time)

## Output result ##
plt.figure()
plt.plot(tim,wave_vel)
plt.plot(tim,output_vel)
plt.show()
