import numpy as np


area_x = 100.0
area_z = 10.0

nx = 20
nz = 5
dof = 2

xg = np.linspace(0,area_x,2*nx+1,endpoint=True)
zg = np.linspace(0,area_z,2*nz+1,endpoint=True)

inode = 0
node = np.empty([len(xg),len(zg)],dtype=np.int32)
node_lines = []
for k in range(len(zg)):
    for i in range(len(xg)):
        dofx,dofz = 1,1
        if k == len(zg)-1:
            dofz = 0
        if i == 0:
            dofz = 0
        if i == len(xg)-1:
            dofz = 0

        node[i,k] = inode
        node_lines += [ "{} {} {} {} {} \n".format(inode,xg[i],zg[k],dofx,dofz) ]
        inode += 1

element_lines = []
ielem = 0
for k in range(nz):
    for i in range(nx):
        style = "2d9solid"

        vs = 250.0
        vp = 1500.0
        rho = 1750.0

        param_line = "{} {} {} {} {} ".format(ielem,style,vs,vp,rho)
        style_line = "{} {} {} {} {} {} {} {} {}".format(node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
                                                         node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
                                                         node[2*i+1,2*k+1])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

for i in range(nx):
    style = "1d3input"

    vs = 250.0
    vp = 1500.0
    rho = 1750.0

    param_line = "{} {} {} {} {} ".format(ielem,style,vs,vp,rho)
    style_line = "{} {} {} ".format(node[2*i,-1],node[2*i+2,-1],node[2*i+1,-1])

    element_lines += [param_line + style_line + "\n"]
    ielem += 1

nnode = inode
nelem = ielem


with open("mesh.in","w") as f:
    f.write("{} {} {} \n".format(nnode,nelem,dof))
    f.writelines(node_lines)
    f.writelines(element_lines)
