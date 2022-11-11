import numpy as np
import matplotlib.pyplot as plt

output_dir = "result/"

d = np.loadtxt(output_dir+'output_element_e.stress_xz').T
tim,d1 = d[0],d[1]
fig,ax = plt.subplots()
ax.set_xlabel('time')
ax.set_ylabel('stress xz', labelpad=4.0)
ax.plot(tim,d1) #label
fig.savefig(output_dir+'xz_e.png')
plt.close(fig)

d = np.loadtxt(output_dir+'result_e.acc').T
d2 = np.loadtxt(output_dir+'input_e.acc').T
tim,d1 = d[0],d[1]
fig,ax = plt.subplots()
ax.set_xlabel('time')
ax.set_ylabel('result acc', labelpad=4.0)
ax.plot(tim,d2[1],label='input acc')
ax.plot(tim,d1,label='output acc') #label
ax.legend()
fig.savefig(output_dir+'acc_e.png')
plt.close(fig)

d = np.loadtxt(output_dir+'result_e.disp').T
tim,d1,d2 = d[0],d[1],d[2]
fig,ax = plt.subplots()
ax.set_xlabel('time')
ax.set_ylabel('result disp', labelpad=4.0)
ax.plot(tim,d1) #label
ax.plot(tim,d2) #label
fig.savefig(output_dir+'disp_e.png')
plt.close(fig)