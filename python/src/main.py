import numpy as np

import io_data
import input_wave
import plot_model

fem = io_data.input_mesh("input/mesh.in")
fem.set_init()
# plot_model.plot_mesh(fem)

fsamp = 5000
tim = np.linspace(0,0.5,int(fsamp*0.5),endpoint=False)
dt = 1.0/fsamp

fem.update_init(dt)
ax = plot_model.plot_mesh_update_init()

vel0 = 0.0
for it in range(len(tim)):
    acc0 = input_wave.ricker(tim[it],0.2,4.0,1.0)
    vel0 += acc0*dt
    fem.update_time(dt,acc0,vel0,input_wave=True)

    if it%50 == 0:
        plot_model.plot_mesh_update(ax,fem,1000.)
        print(it,fem.nodes[20].u)
