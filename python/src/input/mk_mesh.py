import numpy as np
import os

area_x = 100.0
area_z = 10.0

nx = 20
nz = 12
dof = 2

xg = np.linspace(0,area_x,2*nx+1,endpoint=True)
zg = np.linspace(0,area_z,2*nz+1,endpoint=True)


######
width_box = 10.0  # 躯体の幅(m)
nx_box = 4  # 躯体の要素数（水平）
nz_box = 6  # 躯体の要素数（鉛直）

i0_box = (nx-nx_box)//2     # 躯体左端の要素位置
i1_box = (nx+nx_box)//2-1   # 躯体右端の要素位置
k1_box = nz_box-1           # 躯体下端の要素位置

i0_box_node = i0_box * 2         # 躯体左端のノード位置
i1_box_node = (i1_box+1) * 2     # 躯体右端のノード位置
k1_box_node = (k1_box+1) * 2     # 躯体下端のノード位置

#######  log2 で 要素幅を変える
xc = area_x/2.0
x0 = width_box/2.0
x1 = area_x - xc

x0_log2 = np.log(x0)/np.log(2)
x1_log2 = np.log(x1)/np.log(2)
num_log = 2*nx+1 - i1_box_node
log_grid = np.logspace(x0_log2,x1_log2,num_log,base=2)
box_grid = np.linspace(0,x0,nx_box,endpoint=False)

xg_log = np.zeros_like(xg)
xg_log[nx:i1_box_node] = box_grid + xc
xg_log[i1_box_node:]    =  log_grid + xc
xg_log[i0_box_node+1:nx+1] = -box_grid[::-1] + xc
xg_log[0:i0_box_node+1] = -log_grid[::-1] + xc

xg = np.copy(xg_log)
print(xg)

### Set node ###
node = np.empty([len(xg),len(zg)],dtype=np.int32)
node_lines = []

inode = 0
for k in range(len(zg)):
    for i in range(len(xg)):
        dofx,dofz = 1,1
        if k == len(zg)-1:
            dofz = 0

        if i0_box_node <= i <= i1_box_node:
            if k == k1_box_node:
                dofz = 0

        node[i,k] = inode
        node_lines += [ "{} {} {} {} {} \n".format(inode,xg[i],zg[k],dofx,dofz)]
        inode += 1

# 躯体周辺にノードを追加 #
node_double = np.empty([len(xg),len(zg)],dtype=np.int32)

for k in range(k1_box_node+1):
    dofx,dofz = 1,1
    if k == k1_box_node:
        dofz = 0

    i = i0_box_node
    node_double[i,k] = inode
    node_lines += [ "{} {} {} {} {} \n".format(inode,xg[i],zg[k],dofx,dofz)]
    inode += 1

    i = i1_box_node
    node_double[i,k] = inode
    node_lines += [ "{} {} {} {} {} \n".format(inode,xg[i],zg[k],dofx,dofz)]
    inode += 1

for i in range(i0_box_node+1,i1_box_node):
    # dofx,dofz = 1,1
    dofx,dofz = 1,0

    k = k1_box_node
    node_double[i,k] = inode
    node_lines += [ "{} {} {} {} {} \n".format(inode,xg[i],zg[k],dofx,dofz)]
    inode += 1


### Set element ###
element_lines = []

ielem = 0
for k in range(nz):
    for i in range(nx):
        im = 1              # 地盤の材料パラメータ
        if k <= k1_box:
            if i0_box <= i and i <= i1_box:  # 躯体の材料パラメータ
                im = 0

        n0,n1,n2,n3 = node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2]
        n4,n5,n6,n7 = node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1]
        n8 = node[2*i+1,2*k+1]

        if k <= k1_box and i == i0_box:                     # 躯体左端
            n0,n3,n7 = node_double[2*i,2*k],node_double[2*i,2*k+2],node_double[2*i,2*k+1]
        if k <= k1_box and i == i1_box:                     # 躯体右端
            n1,n2,n5 = node_double[2*i+2,2*k],node_double[2*i+2,2*k+2],node_double[2*i+2,2*k+1]
        if k == k1_box and i0_box <= i and i <= i1_box:     # 躯体下端
            n2,n3,n6 = node_double[2*i+2,2*k+2],node_double[2*i,2*k+2],node_double[2*i+1,2*k+2]

        style = "2d9solid"
        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {} {} {} {} {} {}".format(n0,n1,n2,n3,n4,n5,n6,n7,n8)

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

######
for k in range(k1_box_node):          # 側面
    style = "spring"
    im = 3

    param_line = "{} {} {} ".format(ielem,style,im)

    i = i0_box_node
    style_line = "{} {}".format(node[i,k],node_double[i,k])
    element_lines += [param_line + style_line + "\n"]
    ielem += 1

    param_line = "{} {} {} ".format(ielem,style,im)

    i = i1_box_node
    style_line = "{} {}".format(node[i,k],node_double[i,k])
    element_lines += [param_line + style_line + "\n"]
    ielem += 1


for i in range(i0_box_node,i1_box_node+1):         # 下面
    style = "spring"
    im = 4

    param_line = "{} {} {} ".format(ielem,style,im)

    k = k1_box_node
    style_line = "{} {}".format(node[i,k],node_double[i,k])
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
material_lines += ["{} {} {} {} {} \n".format(0,"nu_E_rho",0.33,2*3.904e9*(1+0.33),1700.0*0.5)]
material_lines += ["{} {} {} {} {} \n".format(1,"nu_vs_rho",0.33,150.0,1700.0)]
material_lines += ["{} {} {} {} {} \n".format(2,"nu_vs_rho",0.33,300.0,1700.0)]

material_lines += ["{} {} {} {} {} {} \n".format(3,"spring_normal",1.0e10,1.0e-8,1.0,0.0)] # 側面バネ（法線ベクトル：(1,0)）
material_lines += ["{} {} {} {} {} {} \n".format(4,"spring_normal",1.0e10,1.0e10,0.0,1.0)] # 下面バネ（法線ベクトル：(0,1)）

nmaterial = len(material_lines)

### Set output ###
output_node_lines = []
for i in range(len(xg)):
    output_node_lines += ["{} \n".format(i)]        #define output nodes

output_element_lines = []
ielem = 0
for k in range(nz):
    for i in range(nx):
        if i == i1_box+1:
            output_element_lines += ["{} \n".format(ielem)]
        ielem += 1
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
