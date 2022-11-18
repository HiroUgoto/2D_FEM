import numpy as np
import os,sys

area_x = 1.0
area_z = 10.0

nx = 1
nz = 10
dof = 2

xg = np.linspace(0,area_x,nx+1,endpoint=True)
zg = np.linspace(0,area_z,nz+1,endpoint=True)

### Set node ###
node = np.empty([len(xg),len(zg)],dtype=np.int32)
node_lines = []

inode = 0
for k in range(len(zg)):
    for i in range(len(xg)):
        dofx,dofz = 1,1
        if k == len(zg)-1:
            dofz = 0

        node[i,k] = inode
        node_lines += [ "{} {} {} {} {} \n".format(inode,xg[i],zg[k],dofx,dofz)]
        inode += 1

### Set element ###
element_lines = []

ielem = 0
for k in range(nz):
    for i in range(nx):
        im = 1

        n0,n1,n2,n3 = node[i,k],node[i+1,k],node[i+1,k+1],node[i,k+1]

        style = "2d4solid"
        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(n0,n1,n2,n3)

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

######
for i in range(nx):
    im = 2

    style = "1d2input"
    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {}".format(node[i,-1],node[i+1,-1])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

######
for k in range(len(zg)):     #connected element
    style = "connect"
    im = -1

    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {}".format(node[0,k],node[nx,k])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

nnode = inode       #number of nodes
nelem = ielem       #number of elements


### Set material ###
# dl_style = 'ep_DL1d'          # DeepLearning (full sample length)
dl_style = 'ep_DL1d_light'      # DeepLearning (1000 sample length)
# dl_style = 'ep_DL1d_ghe'      # GHE model (no DeepLearning)

info = {  # soil 10 of cycle data
    'H':17.23,  # m
    'P0':280.0,  # kN/m2
    'N':4.6,  # times
    'G0':48.3e6,  # N/m2
    'sand':3.0,  # %
    'silt':46.0,  # %
    'clay':51.0,  # %
    'wL':47.3,  # negative number instead of np.nan
    'wP':22.0,  # negative number instead of np.nan
}

# alternative elastic parameter
G0e = info['G0']
nu,rho = 0.33,1700
Ee = 2*G0e*(1+nu)

material_lines = []
material_lines += ["{} {} {} {} {} \n".format(0,"nu_E_rho",nu,Ee,rho)]
material_lines += ["{} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(1,dl_style,rho,nu,info['H'],info['P0'],info['N'],info['G0'],info['sand'],info['silt'],info['clay'],info['wL'],info['wP'])]
material_lines += ["{} {} {} {} {}\n".format(2,"nu_vs_rho",0.33,300.0,1700.0)]

nmaterial = len(material_lines)

### Set output ###
output_node_lines = []
for i in range(len(xg)):
    output_node_lines += ["{}\n".format(i)]

output_element_lines = []
ielem = 0
for k in range(nz):
    for i in range(nx):
        output_element_lines += ["{}\n".format(ielem)]
        ielem += 1

output_nnode = len(output_node_lines)
output_nelem = len(output_element_lines)

with open('mesh.in',"w") as f:
    f.write("{} {} {} {}\n".format(nnode,nelem,nmaterial,dof))
    f.writelines(node_lines)
    f.writelines(element_lines)
    f.writelines(material_lines)

with open('output.in',"w") as f:
    f.write("{} {}\n".format(output_nnode,output_nelem))
    f.writelines(output_node_lines)
    f.writelines(output_element_lines)
