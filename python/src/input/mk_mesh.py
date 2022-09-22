import numpy as np
import os

area_x = 6.0
area_z = 6.0

nx = 3
nz = 3
dof = 2

xg = np.linspace(0,area_x,nx+1,endpoint=True)
zg = np.linspace(0,area_z,nz+1,endpoint=True)


### Set node ###
node = np.empty([len(xg),len(zg)],dtype=np.int32)       #node_idを振った配列(転置)
node_lines = []

inode = 0
for k in range(len(zg)):
    for i in range(len(xg)):
        dofx,dofz = 1,1
        if i == 0:
            dofx = 0
        elif i == len(xg)-1:
            dofx = 0
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
        im = 0
        if i == 1:
            style = "2d4solidX"
        else:
            style = "2d4solid"

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(node[i,k],node[i+1,k],node[i+1,k+1],node[i,k+1])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

nnode = inode       #number of nodes
nelem = ielem       #number of elements


### Set material ###
material_lines = []
material_lines += ["{} {} {} {} {} \n".format(0,"vs_vp_rho",100.0,1500.0,1750.0)]

nmaterial = len(material_lines)

### Set output ###
output_node_lines = []
for i in range(len(xg)):
    output_node_lines += ["{} \n".format(i)]        #define output nodes

output_element_lines = []
# for i in range(0,nelem-nx-len(zg)):        #define output elements
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

# with open("var.in","w") as f:       #save var, depend on target area
#     f.write("{} {}\n".format(modelid,"modelid"))
#     f.write("{} {} {} {}\n".format(area_x,area_z,"area_x","area_z"))
#     f.write("{} {} {} {}\n".format(nx,nz,"nx","nz"))
#     if modelid == 1:
#         f.write("{} {} {} {}\n".format(nx1,nx2,"nx1","nx2"))
#         f.write("{} {} {} {}\n".format(nz1,nz2,"nz1","nz2"))
