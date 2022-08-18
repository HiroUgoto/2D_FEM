import numpy as np
import os

area_x = 5000.0
area_z = 5000.0

nx = 10
nz = 10
dof = 2

xg = np.linspace(0,area_x,2*nx+1,endpoint=True)
zg = np.linspace(0,area_z,2*nz+1,endpoint=True)


### Set node ###
node = np.empty([len(xg),len(zg)],dtype=np.int32)       #node_idを振った配列(転置)
node_lines = []

inode = 0
for k in range(len(zg)):
    for i in range(len(xg)):
        dofx,dofz = 1,1
        # if k == len(zg)-1:
        #     dofz = 0

        node[i,k] = inode
        node_lines += [ "{} {} {} {} {} \n".format(inode,xg[i],zg[k],dofx,dofz)]
        inode += 1


### Set element ###
element_lines = []

ielem = 0
for k in range(nz):
    for i in range(nx):
        im = 0

        style = "2d9solid"

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {} {} {} {} {} {}".format(node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
                                                         node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
                                                         node[2*i+1,2*k+1])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

for k in range(nz):
    im = 0
    style = "1d3visco"

    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {} {}".format(node[0,2*k],node[0,2*k+2],node[0,2*k+1])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {} {}".format(node[2*nx,2*k+2],node[2*nx,2*k],node[2*nx,2*k+1])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

for i in range(nx):
    im = 0
    style = "1d3visco"

    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {} {}".format(node[2*i,2*nz],node[2*i+2,2*nz],node[2*i+1,2*nz])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

nnode = inode       #number of nodes
nelem = ielem       #number of elements


### Set material ###
material_lines = []
material_lines += ["{} {} {} {} {} \n".format(0,"vs_vp_rho",3464.0,6000.0,2450.0)]

nmaterial = len(material_lines)

### Set output ###
output_node_lines = []
output_node_lines += ["{} \n".format(node[nx,0])]
# for i in range(len(xg)):
#     output_node_lines += ["{} \n".format(i)]

output_element_lines = []
# for i in range(0,nelem-nx-len(zg)):
#     output_element_lines += ["{} \n".format(i)]

output_nnode = len(output_node_lines)
output_nelem = len(output_element_lines)


with open("mesh.in","w") as f:
    f.write("{} {} {} {} \n".format(nnode,nelem,nmaterial,dof))
    f.writelines(node_lines)
    f.writelines(element_lines)
    f.writelines(material_lines)

with open("output.in","w") as f:
    f.write("{} {} \n".format(output_nnode,output_nelem))
    f.writelines(output_node_lines)
    f.writelines(output_element_lines)
