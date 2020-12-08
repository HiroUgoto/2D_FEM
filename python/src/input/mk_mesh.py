import numpy as np

modelid = 1     #0:square mesh,1:flexible mesh

### Set target area ###
if  modelid == 0:
    area_x = 50.0
    area_z = 10.0

    nx =  20
    nz =  5
    dof = 2

    xg = np.linspace(0,area_x,2*nx+1,endpoint=True)
    zg = np.linspace(0,area_z,2*nz+1,endpoint=True)

elif modelid == 1:
    area_x = 70.0
    area_z = 5.0

    nx1,nx2 = 4, 8
    nz1,nz2 = 3, 2
    nx =  nx1 + nx2 + nx1
    nz =  nz1 + nz2
    dof = 2

    zg = np.linspace(0,area_z,2*nz+1,endpoint=True)
    xg = np.empty([2*nx+1,2*nz+1])       #全nodex座標

    for k in range(2*nz1+1):
        xg[:2*nx1,k] = np.linspace(0,20+10/(len(zg)-1)*k,2*nx1,endpoint=False)
        xg[2*nx1:2*(nx1+nx2),k] = np.linspace(20+10/(len(zg)-1)*k,area_x-(20+10/(len(zg)-1)*k),2*nx2,endpoint=False)
        xg[2*(nx1+nx2):,k] = np.linspace(area_x-(20+10/(len(zg)-1)*k),area_x,2*nx1+1,endpoint=True)      #傾斜部分

    for k in range(2*nz1+1,2*nz+1):
        xg[:,k] = np.copy(xg[:,2*nz1])


### Set node ###
inode = 0

if modelid == 0:
    node = np.empty([len(xg),len(zg)],dtype=np.int32)       #node_idを振った配列(転置)
    node_lines = []
    for k in range(len(zg)):
        for i in range(len(xg)):
            dofx,dofz = 1,1
            dofx_static,dofz_static = 1,1
            if k == len(zg)-1:
                dofz = 0
                dofz_static = 0
            if i == 0:
                dofz = 0
                dofx_static = 0
            if i == len(xg)-1:
                dofz = 0
                dofx_static = 0

            node[i,k] = inode
            node_lines += [ "{} {} {} {} {} {} {}\n".format(inode,xg[i],zg[k],dofx,dofz,dofx_static,dofz_static) ]
            inode += 1

elif modelid == 1:
    node = np.empty([len(xg[:,0]),len(zg)],dtype=np.int32)       #node_idを振った配列(転置)
    node_lines = []
    for k in range(len(zg)):
        for i in range(len(xg[:,0])):
            dofx,dofz = 1,1
            if k == len(zg)-1:
                dofz = 0

            node[i,k] = inode
            node_lines += [ "{} {} {} {} {}\n".format(inode,xg[i,k],zg[k],dofx,dofz) ]
            inode += 1

### Set element ###
element_lines = []
ielem = 0

if modelid == 0:
    for k in range(nz):
        im = 1
        for i in range(nx):
            im = 1
            if k <= 1 and i >= int(nx/2):
                im = 0

            style = "2d9solid"

            param_line = "{} {} {} ".format(ielem,style,im)
            style_line = "{} {} {} {} {} {} {} {} {}".format(node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
                                                             node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
                                                             node[2*i+1,2*k+1])

            element_lines += [param_line + style_line + "\n"]
            ielem += 1

elif modelid == 1:
    for k in range(nz):
        im = 1
        for i in range(nx):
            im = 1
            if k < nz1 and nx1 <= i < nx1+nx2:
                im = 0

            style = "2d9solid"

            param_line = "{} {} {} ".format(ielem,style,im)
            style_line = "{} {} {} {} {} {} {} {} {}".format(node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
                                                             node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
                                                             node[2*i+1,2*k+1])

            element_lines += [param_line + style_line + "\n"]
            ielem += 1


for i in range(nx):
    style = "1d3input"
    im = 1

    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {} {} ".format(node[2*i,-1],node[2*i+2,-1],node[2*i+1,-1])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1


for k in range(len(zg)):     #connected element
    style = "connect"
    im = -1

    param_line = "{} {} {} ".format(ielem,style,im)
    style_line = "{} {}".format(node[0,k],node[2*nx,k])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1


nnode = inode       #nodeの総数
nelem = ielem       #elementの総数


### Set material ###
material_lines = []
material_lines += ["{} {} {} {} {} \n".format(0,"vs_vp_rho",0.0,1500.0,1750.0)]
material_lines += ["{} {} {} {} {} \n".format(1,"vs_vp_rho",200.0,1500.0,1750.0)]

nmaterial = len(material_lines)


### Set output ###
output_node_lines = []
for i in range(0,nnode):
    output_node_lines += ["{} \n".format(i)]       #outputするnode指定

output_element_lines = []
for i in range(0,nelem-nx-len(zg)):        #outputするelement番号指定
    output_element_lines += ["{} \n".format(i)]

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
