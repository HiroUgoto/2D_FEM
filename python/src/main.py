import numpy as np

import io_data
import input_wave
import plot_model

fem = io_data.input_mesh("input/mesh.in")
fem.set_init()
# plot_model.plot_mesh(fem)


fsamp = 2500
tim = np.linspace(0,1.0,fsamp,endpoint=False)
dt = 1.0/fsamp

ax = plot_model.plot_mesh_update_init()

for it in range(len(tim)):
    acc0 = input_wave.ricker(tim[it],0.25,4.0,1.0)
    fem.update_time(dt,acc0)

    if it%25 == 0:
        plot_model.plot_mesh_update(ax,fem,1000.)
        print(it,fem.nodes[20].u)
