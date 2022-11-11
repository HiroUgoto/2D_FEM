import numpy as np
import os,sys

area_x = 1.0
area_z = 10.0

nx = 1
nz = 10
dof = 2

xg = np.linspace(0,area_x,2*nx+1,endpoint=True)
zg = np.linspace(0,area_z,2*nz+1,endpoint=True)

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

        n0,n1,n2,n3 = node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2]
        n4,n5,n6,n7 = node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1]
        n8 = node[2*i+1,2*k+1]

        style = "2d9solid"
        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {} {} {} {} {} {}".format(n0,n1,n2,n3,n4,n5,n6,n7,n8)

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

######
for i in range(nx):
    style = "1d3input"
    im = 2

    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {} {} ".format(node[2*i,-1],node[2*i+2,-1],node[2*i+1,-1])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

######
for k in range(len(zg)):     #connected element
    style = "connect"
    im = -1

    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {}".format(node[0,k],node[2*nx,k])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

nnode = inode       #number of nodes
nelem = ielem       #number of elements


### Set material ###
material_lines = []
material_lines += ["{} {} {} {} {} \n".format(0,"nu_E_rho",0.2,100.0e9,850.0)]
# material_lines += ["{} {} {} {} {} \n".format(1,"nu_vs_rho",0.33,150.0,1700.0)]
# material_lines += ["{} {} {} {} {} {} {} {} {} {}\n".format(1,"ep_Li",1700.0,0.33,420,0.97,0.7148,0.957,0.0,4.e3)]
#                                                      # rho, nu, G0, M, e0, eg, d1, cohesion
material_lines += ["{} {} {} {} {} {} {} {} {} {}\n".format(1,"ep_eff_Li",1700.0,0.33,420,0.97,0.7148,0.957,0.41,4.e3)]
                                                     # rho, nu, G0, M, e0, eg, d1, cohesion
material_lines += ["{} {} {} {} {} \n".format(2,"nu_vs_rho",0.33,350.0,1800.0)]


nmaterial = len(material_lines)

### Set output ###
output_node_lines = []
for i in range(len(xg)):
    output_node_lines += ["{} \n".format(i)]        #define output nodes

output_element_lines = []
ielem = 0
for k in range(nz):
    for i in range(nx):
        output_element_lines += ["{} \n".format(ielem)]
        ielem += 1

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
