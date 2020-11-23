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
output_velx = np.zeros((ntim,fem.output_nnode))
output_velz = np.zeros((ntim,fem.output_nnode))
output_strainxx = np.zeros((ntim,fem.output_nelem))
output_strainyy = np.zeros((ntim,fem.output_nelem))
output_strainxy = np.zeros((ntim,fem.output_nelem))
output_dispx = np.zeros((ntim,fem.output_nnode))
output_dispz = np.zeros((ntim,fem.output_nnode))

for it in range(len(tim)):
    acc0 = np.array([0.0,0.0])
    vel0 = np.array([wave_vel[it],0.0])

    fem.update_time(acc0,vel0,input_wave=True)

    output_velx[it,:] = [node.v[0] for node in fem.output_nodes]     #axis [0]:x [1]:y
    output_velz[it,:] = [node.v[1] for node in fem.output_nodes]
    output_strainxx[it,:] = [element.strain[0] for element in fem.output_elements]        #axis [0]:xx [1]:yy [2]:xy
    output_strainyy[it,:] = [element.strain[1] for element in fem.output_elements]        #axis [0]:xx [1]:yy [2]:xy
    output_strainxy[it,:] = [element.strain[2] for element in fem.output_elements]        #axis [0]:xx [1]:yy [2]:xy
    output_dispx[it,:] = [node.u[0] for node in fem.output_nodes]
    output_dispz[it,:] = [node.u[1] for node in fem.output_nodes]

    if it%50 == 0:      #terminal出力の時間間隔
        plot_model.plot_mesh_update(ax,fem,500.)
        plt.savefig("output_fig\img_"+str(it)+".png")
        print(it,output_velx[it,0],output_strainxx[it,0])

## --- Write output file --- ##
output_tim = np.arange(ntim).reshape(ntim,1)

output_w_velx = np.hstack((output_tim,output_velx))
np.savetxt("output\output_velx.dat",output_w_velx,delimiter="    ")
output_w_velz = np.hstack((output_tim,output_velz))
np.savetxt("output\output_velz.dat",output_w_velz,delimiter="    ")

output_w_strainxx = np.hstack((output_tim,output_strainxx))
np.savetxt("output\output_strainxx.dat",output_w_strainxx,delimiter="    ")
output_w_strainyy = np.hstack((output_tim,output_strainyy))
np.savetxt("output\output_strainyy.dat",output_w_strainyy,delimiter="    ")
output_w_strainxy = np.hstack((output_tim,output_strainxy))
np.savetxt("output\output_strainxy.dat",output_w_strainxy,delimiter="    ")

output_w_dispx = np.hstack((output_tim,output_dispx))
np.savetxt("output\output_dispx.dat",output_w_dispx,delimiter="    ")
output_w_dispz = np.hstack((output_tim,output_dispz))
np.savetxt("output\output_dispz.dat",output_w_dispz,delimiter="    ")

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

## Output result ##
plt.figure()
plt.plot(tim,wave_vel)
plt.plot(tim,output_velx[:,0])
plt.show()
