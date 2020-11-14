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
plot_model.plot_mesh(fem)


## --- Define input wave --- ##
fsamp = 5000
duration = 0.5

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
wave_acc = input_wave.tapered_sin(tim,fp=5.0,taper=0.2,duration=1.0,amp=1.0)
wave_vel = np.cumsum(wave_acc) * dt
ntim = len(tim)

## --- Prepare time solver --- ##
fem.update_init(dt)
ax = plot_model.plot_mesh_update_init()

## Iteration ##
output_vel = np.zeros((ntim,fem.output_nnode))
output_strain = np.zeros((ntim,fem.output_nelem))
output_disp = np.zeros((ntim,fem.output_nnode))

for it in range(len(tim)):
    acc0 = np.array([0.0,0.0])
    vel0 = np.array([wave_vel[it],0.0])

    fem.update_time(acc0,vel0,input_wave=True)

    output_vel[it,:] = [node.v[0] for node in fem.output_nodes]     #axis [0]:x [1]:y
    output_strain[it,:] = [element.strain[0] for element in fem.output_elements]        #axis [0]:xx [1]:yy [2]:xy
    output_disp[it,:] = [node.u[0] for node in fem.output_nodes]


    if it%10 == 0:      #terminal出力の時間間隔
        plot_model.plot_mesh_update(ax,fem,500.)
        print(it,output_vel[it,0],output_strain[it,0])

## --- Write output file --- ##
output_tim = np.arange(ntim).reshape(ntim,1)

output_w_vel = np.hstack((output_tim,output_vel))
np.savetxt("output\output_vel.dat",output_w_vel,delimiter="    ")
output_w_strain = np.hstack((output_tim,output_strain))
np.savetxt("output\output_strain.dat",output_w_strain,delimiter="    ")
output_w_disp = np.hstack((output_tim,output_disp))
np.savetxt("output\output_disp.dat",output_w_disp,delimiter="    ")

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

## Output result ##
plt.figure()
plt.plot(tim,wave_vel)
plt.plot(tim,output_vel[:,0])
plt.show()
