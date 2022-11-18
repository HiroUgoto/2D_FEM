import numpy as np
import os,sys

plastic = True
# plastic,ep = False,None
ep = 'ep_DL1d'
# ep = 'ep_DL1d_light'
# ep = 'ep_DL1d_ghe'

# area_x = 10.0
# area_z = 10.0
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
            # dofx = 0

        node[i,k] = inode
        node_lines += [ "{} {} {} {} {} \n".format(inode,xg[i],zg[k],dofx,dofz)]
        inode += 1

### Set element ###
element_lines = []

ielem = 0
for k in range(nz):
    for i in range(nx):
        if plastic:
            im = 1              # 地盤の材料パラメータ
        else:  # elastic
            im = 0

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
G0 = info['G0']
nu,rho = 0.33,1700
E = 2*G0*(1+nu)

material_lines = []
# material_lines += ["{} {} {} {} {}\n".format(0,"nu_E_rho",0.33,2*3.904e9*(1+0.33),1700.0*0.5)]
material_lines += [f"0 nu_E_rho {nu} {E} {rho}\n"]
# material_lines += ["{} {} {} {} {} {} {} {} {} {}\n".format(1,"ep_Li",1700.0,0.33,210,0.93,0.7148,0.957,0.0,4.e3)]
# rho, nu, G0, M, e0, eg, d1, cohesion

if plastic:
    material_lines += [f"1 {ep} {rho} {nu} {info['H']} {info['P0']} {info['N']} {info['G0']} {info['sand']} {info['silt']} {info['clay']} {info['wL']} {info['wP']}\n"]
    # rho, nu, info
    material_lines += ["{} {} {} {} {}\n".format(2,"nu_vs_rho",0.33,300.0,1700.0)]


nmaterial = len(material_lines)

### Set output ###
output_node_lines = []
for i in range(len(xg)):
    output_node_lines += ["{}\n".format(i)]        #define output nodes

output_element_lines = []
ielem = 0
for k in range(nz):
    for i in range(nx):
        output_element_lines += ["{}\n".format(ielem)]
        ielem += 1

output_nnode = len(output_node_lines)
output_nelem = len(output_element_lines)


file = 'mesh.in' if plastic else 'mesh_elastic.in'
with open(file,"w") as f:
    f.write("{} {} {} {}\n".format(nnode,nelem,nmaterial,dof))
    f.writelines(node_lines)
    f.writelines(element_lines)
    f.writelines(material_lines)

file = 'output.in' if plastic else 'output_elastic.in'
with open(file,"w") as f:
    f.write("{} {}\n".format(output_nnode,output_nelem))
    f.writelines(output_node_lines)
    f.writelines(output_element_lines)
