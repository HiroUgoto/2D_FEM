import matplotlib.pyplot as plt
import numpy as np
import time

import io_data
import input_wave
import plot_model
import os
import shutil
import datetime

dir = "output/outputdir"
os.makedirs(dir,exist_ok=True)      #make output folder
os.makedirs(dir+"/fig",exist_ok=True)

start = time.time()
## --- Input FEM Mesh --- ##
fem = io_data.input_mesh("input/mesh.in")
outputs = io_data.input_outputs("input/output.in")

## --- FEM Set up --- ##
fem.set_init()
fem.set_output(outputs)
# plot_model.plot_mesh(fem)

## --- Define input wave --- ##
fsamp = 1600
duration = 0.5

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
wave_acc = input_wave.ricker(tim,fp=1.0,tp=1.0,amp=2.0)
wave_vel = np.cumsum(wave_acc) * dt
ntim = len(tim)

# plt.figure()
# plt.plot(tim,wave_vel)
# plt.show()

ax = plot_model.plot_mesh_update_init()

## --- Static deformation --- ##
fem.self_gravity()
# plot_model.plot_mesh_update(ax,fem,10.)

## --- Prepare time solver --- ##
fem.update_init(dt)

## Iteration ##
output_velx = np.zeros((ntim,fem.output_nnode))
output_velz = np.zeros((ntim,fem.output_nnode))
output_strainxx = np.zeros((ntim,fem.output_nelem))
output_strainzz = np.zeros((ntim,fem.output_nelem))
output_strainxz = np.zeros((ntim,fem.output_nelem))
output_dispx = np.zeros((ntim,fem.output_nnode))
output_dispz = np.zeros((ntim,fem.output_nnode))

for it in range(len(tim)):
    acc0 = np.array([0.0,0.0])
    vel0 = np.array([wave_vel[it],0.0])

    fem.update_time(acc0,vel0,input_wave=True,FD=True)
    # fem.update_time(acc0,vel0,input_wave=True)

    output_dispx[it,:] = [node.u[0]-node.u0[0] for node in fem.output_nodes]
    output_dispz[it,:] = [node.u[1]-node.u0[1] for node in fem.output_nodes]
    output_velx[it,:] = [node.v[0] for node in fem.output_nodes]     #axis [0]:x [1]:y
    output_velz[it,:] = [node.v[1] for node in fem.output_nodes]
    output_strainxx[it,:] = [element.strain[0] for element in fem.output_elements]        #axis [0]:xx [1]:zz [2]:xz
    output_strainzz[it,:] = [element.strain[1] for element in fem.output_elements]
    output_strainxz[it,:] = [element.strain[2] for element in fem.output_elements]

    if it%20 == 0:
        # plot_model.plot_mesh_update(ax,fem,10.)
        # print(it,"t=",it*dt,output_dispx[it,:])
        print(it,"t=",it*dt,output_dispx[it,8])
        plt.savefig(dir+"/fig/img_"+str(it)+".png")


# plot_model.plot_mesh_update(ax,fem,10.,fin=True)

## --- Write output file --- ##
# with open("input/var.in","a") as f:
#     f.write("{} {} {} {}\n".format("inputwave",fsamp,duration,)
# shutil.copy("input/mesh.in",dir)        #movefile to output folder
# shutil.copy("input/output.in",dir)
# shutil.copy("input/var.in",dir)
#
# output_tim = np.arange(ntim).reshape(ntim,1)
#
# output_w_velx = np.hstack((output_tim,output_velx))
# np.savetxt(dir+"/velx.dat",output_w_velx,delimiter="    ")
# output_w_velz = np.hstack((output_tim,output_velz))
# np.savetxt(dir+"/velz.dat",output_w_velz,delimiter="    ")
#
# output_w_strainxx = np.hstack((output_tim,output_strainxx))
# np.savetxt(dir+"/strainxx.dat",output_w_strainxx,delimiter="    ")
# output_w_strainzz = np.hstack((output_tim,output_strainzz))
# np.savetxt(dir+"/strainzz.dat",output_w_strainzz,delimiter="    ")
# output_w_strainxz = np.hstack((output_tim,output_strainxz))
# np.savetxt(dir+"/strainxz.dat",output_w_strainxz,delimiter="    ")
#
# output_w_dispx = np.hstack((output_tim,output_dispx))
# np.savetxt(dir+"/dispx.dat",output_w_dispx,delimiter="    ")
# output_w_dispz = np.hstack((output_tim,output_dispz))
# np.savetxt(dir+"/dispz.dat",output_w_dispz,delimiter="    ")
#
# now = datetime.datetime.now()
# os.rename(dir,"output/"+now.strftime("%Y%m%d-%H%M"))        #"outputdir"-->>"time"

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# ## Output result ##
# plt.figure()
# plt.plot(tim,wave_vel)
# plt.plot(tim,output_velx[:,0])
# plt.plot(tim,output_velx[:,6])
# plt.show()
