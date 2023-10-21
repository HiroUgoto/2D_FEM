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
output_dir = "result/"
datasamplerate = 2000  #sampling rate[Hz]
# datasamplerate = 100  #sampling rate[Hz]

## --- FEM Set up --- ##
fem.set_init()
fem.set_output(outputs)
# plot_model.plot_mesh(fem)

## --- Define input wave --- ##
fp = 3.3
amp = 1.e-5
print("Input frequency(Hz):",fp,"Input amplitude(m/s2):",amp)

## --- EP Set up --- ##
fem.set_ep_initial_state()
# fem.set_rayleigh_damping(fp,10*fp,0.002)

## --- Define input wave --- ##
fsamp = 2000
duration = 3.0/fp + 1.0/fp

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
# wave_acc = input_wave.tapered_sin(tim,fp,3.0/fp,duration-1.0/fp,amp)
forced_disp = input_wave.tapered_sin(tim,fp,3.0/fp,duration-1.0/fp,amp) #FEM試験用
forced_nodes = [2,5,8] #FEM試験用(強制変位S)
ntim = len(tim)

# plt.figure()
# plt.plot(tim,wave_acc)
# plt.show()

## --- Prepare time solver --- ##
# ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Time iteration ##
def writef(line,outf,tim,data):
    line[0] = tim
    line[1:] = data
    np.savetxt(outf,line[np.newaxis])
    outf.flush()

def writef_in(outf,data):
    np.savetxt(outf,data)
    outf.flush()

acc0 = np.zeros([2])
vel0 = np.zeros([2])
dis0 = np.zeros([2])

output_accx = np.zeros(fem.output_nnode+1)
output_accz = np.zeros(fem.output_nnode+1)
output_velx = np.zeros(fem.output_nnode+1)
output_velz = np.zeros(fem.output_nnode+1)
output_dispx = np.zeros(fem.output_nnode+1)
output_dispz = np.zeros(fem.output_nnode+1)
output_element_stress_xx = np.zeros(fem.output_nelem+1)
output_element_stress_zz = np.zeros(fem.output_nelem+1)
output_element_stress_xz = np.zeros(fem.output_nelem+1)
output_element_stress_yy = np.zeros(fem.output_nelem+1)
output_element_strain_xx = np.zeros(fem.output_nelem+1)
output_element_strain_zz = np.zeros(fem.output_nelem+1)
output_element_strain_xz = np.zeros(fem.output_nelem+1)
output_element_pw = np.zeros(fem.output_nelem+1)

# accx = open(output_dir+"x.acc","w") 
# accz = open(output_dir+"z.acc","w") 
# velx = open(output_dir+"x.vel","w") 
# velz = open(output_dir+"z.vel","w") 
# dispx = open(output_dir+"x.disp","w") 
# dispz = open(output_dir+"z.disp","w") 
# strainxx = open(output_dir+"xx.str","w")
# strainzz = open(output_dir+"zz.str","w")
# strainxz = open(output_dir+"xz.str","w")
# stressxx = open(output_dir+"stressxx","w")
# stresszz = open(output_dir+"stresszz","w")
# stressxz = open(output_dir+"stressxz","w")
# stressyy = open(output_dir+"stressyy","w")
# stresspw = open(output_dir+"stresspw","w")
# accin = open(output_dir+"acc.in","w") 
# velin = open(output_dir+"vel.in","w") 
# dispin = open(output_dir+"disp.in","w") 

# for it in range(ntim):
for it in range(25):
    # acc0 = np.array([wave_acc[it],0.0])
    # vel0 += acc0*dt
    # dis0 += vel0*dt
    forced_disp0 = np.array([forced_disp[it],0.0]) #FEM試験用

    # fem.update_time(acc0)
    # fem.update_time(acc0,FD=True)
    # fem.update_time(acc0,vel0,input_wave=True)
    fem.update_time_disp(forced_disp0,forced_nodes)


    if it%(fsamp/datasamplerate) == 0:
        # writef_in(accin, [acc0[0]])
        # writef_in(velin, [vel0[0]])
        # writef_in(dispin, [dis0[0]])
        # writef(output_accx, accx, tim[it], [node.a[0] for node in fem.output_nodes] + acc0[0])
        # writef(output_accz, accz, tim[it], [node.a[1] for node in fem.output_nodes])
        # writef(output_dispx, dispx, tim[it], [node.u[0] for node in fem.output_nodes])
        # writef(output_dispz, dispz, tim[it], [node.u[1] for node in fem.output_nodes])
        # writef(output_velx, velx, tim[it], [node.v[0] for node in fem.output_nodes])
        # writef(output_velz, velz, tim[it], [node.v[1] for node in fem.output_nodes])
        # writef(output_element_stress_xx, stressxx, tim[it], [element.stress[0] for element in fem.output_elements])
        # writef(output_element_stress_zz, stresszz, tim[it], [element.stress[1] for element in fem.output_elements])
        # writef(output_element_stress_xz, stressxz, tim[it], [element.stress[2] for element in fem.output_elements])
        # # writef(output_element_stress_xx, stressxx, tim[it], [element.eff_stress[0] for element in fem.output_elements])
        # # writef(output_element_stress_zz, stresszz, tim[it], [element.eff_stress[1] for element in fem.output_elements])
        # # writef(output_element_stress_xz, stressxz, tim[it], [element.eff_stress[2] for element in fem.output_elements])
        # writef(output_element_strain_xx, strainxx, tim[it], [element.strain[0] for element in fem.output_elements])
        # writef(output_element_strain_zz, strainzz, tim[it], [element.strain[1] for element in fem.output_elements])
        # writef(output_element_strain_xz, strainxz, tim[it], [element.strain[2] for element in fem.output_elements])
        # writef(output_element_stress_yy, stressyy, tim[it], [element.stress_yy for element in fem.output_elements])
        # # writef(output_element_stress_yy, stressyy, tim[it], [element.eff_stress_yy for element in fem.output_elements])
        # # writef(output_element_pw, stresspw, tim[it], [element.excess_pore_pressure for element in fem.output_elements])

        print(it,"t=",it*dt,forced_disp0[0],fem.elements[0].stress[0],fem.elements[0].stress[1],fem.elements[0].stress[2],fem.elements[0].stress_yy)  
        print("----------------------------------------------------------------------")


    # if it%(fsamp/datasamplerate*10) == 0:
        # plot_model.plot_mesh_update(ax,fem,100.)    
        # print("t=",it*dt,output_element_stress_zz[1],output_element_pw[1])
        # print(it,"t=",it*dt,output_accx[0],output_element_stress_xx[0],output_element_stress_zz[0],output_element_stress_yy[0])  
        # break

# accin.close()
# velin.close()
# dispin.close()
# accx.close()
# accz.close()
# velx.close()
# velz.close()
# dispx.close()
# dispz.close()
# strainxx.close()
# strainzz.close()
# strainxz.close()
# stressxx.close()
# stresszz.close()
# stressxz.close()
# stressyy.close()
# stresspw.close()

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")
# plot_model.plot_mesh_update(ax,fem,100.,fin=True)
