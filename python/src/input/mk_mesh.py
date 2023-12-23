import numpy as np
import os

area_x = 5000.0
area_z = 3000.0

nx = 200
nz = 120
# dof = 2
dof = 1

# xg = np.linspace(0,area_x,2*nx+1,endpoint=True)
# zg = np.linspace(0,area_z,2*nz+1,endpoint=True)
xg = np.linspace(-1000,area_x-1000,nx+1,endpoint=True)
zg = np.linspace(0,area_z,nz+1,endpoint=True)

### Set node ###
node = np.empty([len(xg),len(zg)],dtype=np.int32) 
node_lines = []

inode = 0
for k in range(len(zg)):
    for i in range(len(xg)):
        # dofx,dofz = 1,1
        dofy = 1

        node[i,k] = inode
        # node_lines += [ "{} {} {} {} {} \n".format(inode,xg[i],zg[k],dofx,dofz)]
        node_lines += [ "{} {} {} {} \n".format(inode,xg[i],zg[k],dofy)]
        inode += 1


### Set element ###
element_lines = []

ielem = 0
for k in range(nz):
    for i in range(nx):
        im = 0

        # style = "2d9solid"
        # param_line = "{} {} {} ".format(ielem,style,im)
        # style_line = "{} {} {} {} {} {} {} {} {}".format(node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
        #                                                  node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
        #                                                  node[2*i+1,2*k+1])

        style = "2d4solid"
        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(node[i,k],node[i+1,k],node[i+1,k+1],node[i,k+1])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

for k in range(nz):
    im = 0
    # style = "1d3visco"
    style = "1d2visco"

    # param_line = "{} {} {} ".format(ielem,style,im)
    # style_line = "{} {} {}".format(node[0,2*k],node[0,2*k+2],node[0,2*k+1])
    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {}".format(node[0,k],node[0,k+1])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

    # param_line = "{} {} {} ".format(ielem,style,im)
    # style_line = "{} {} {}".format(node[2*nx,2*k+2],node[2*nx,2*k],node[2*nx,2*k+1])
    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {}".format(node[nx,k],node[nx,k+1])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

for i in range(nx):
    im = 0
    # style = "1d3visco"
    style = "1d2visco"

    # param_line = "{} {} {} ".format(ielem,style,im)
    # style_line = "{} {} {}".format(node[2*i,2*nz],node[2*i+2,2*nz],node[2*i+1,2*nz])
    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {}".format(node[i,nz],node[i+1,nz])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

nnode = inode       #number of nodes
nelem = ielem       #number of elements


### Set material ###
material_lines = []
material_lines += ["{} {} {} {} {} \n".format(0,"vs_vp_rho",400.0,1500.0,1800.0)]

nmaterial = len(material_lines)

### Set output ###
output_node_lines = []
for i in range(len(xg)):
    output_node_lines += ["{} \n".format(node[i,0])]

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
