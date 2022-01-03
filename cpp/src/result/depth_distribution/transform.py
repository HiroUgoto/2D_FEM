import numpy as np
import matplotlib.pyplot as plt

# nelem = 40
nelem = 2
nnode = nelem

dz = 10.0/nelem

K0 = 0.5
rho = 1700

stress_list_file = "../output_element_list.dat"
stress_file = "../output_element.stress_xx"
acc_file = "../result.acc"

stress_list = np.loadtxt(stress_list_file)
stress_pack = np.loadtxt(stress_file)
acc = np.loadtxt(acc_file,usecols=(1,))

z = stress_list[:,1]
area = np.ones(nnode)*dz
F0 = rho*9.8*z*K0

tim = stress_pack[:,0]
ntim = len(tim)
stress_node = np.zeros([ntim,nnode])
stress = np.zeros(ntim)

for i in range(ntim):
    stress_node[i,:] = -stress_pack[i,1:] + F0[:]

    for inode in range(int(nnode//2)):
        stress[i] += stress_node[i,inode]*area[inode]

output_line = np.vstack([tim,stress]).T
np.savetxt("earth_pressure_time.dat",output_line)


t_act_ind = np.argmin(stress)
t_pas_ind = np.argmax(stress)
t_act = tim[t_act_ind]
t_pas = tim[t_pas_ind]

print(t_act,t_pas)

F_act = stress_node[t_act_ind,:]
F_pas = stress_node[t_pas_ind,:]

acc_act = acc[t_act_ind]
acc_pas = acc[t_pas_ind]
print(acc_act,acc_pas)

output_line = np.vstack([z,F0,F_act,F_pas]).T
np.savetxt("earth_pressure_depth.dat",output_line)

plt.figure()
plt.ylim([10,0])
plt.plot(F0,z,color='k',marker="o")
plt.plot(F_act,z,color='r',marker="o")
plt.plot(F_pas,z,color='b',marker="o")
plt.show()
