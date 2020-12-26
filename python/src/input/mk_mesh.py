import numpy as np
import os

modelid = 2     #0:square mesh,1:flexible mesh

### Set target area ###
## add model---fix make var.in,mk_vtk.py ##
if  modelid == 0:
    area_x = 50.0
    area_z = 10.0

    nx =  20
    nz =  5
    dof = 2

    xg = np.linspace(0,area_x,2*nx+1,endpoint=True)
    zg = np.linspace(0,area_z,2*nz+1,endpoint=True)

elif modelid == 1:
    area_x = 50.0
    area_z = 5.0

    nx1,nx2 = 1, 4
    nz1,nz2 = 1, 1
    nx =  nx1 + nx2 + nx1
    nz =  nz1 + nz2
    dof = 2

    zg = np.linspace(0,area_z,2*nz+1,endpoint=True)
    xg = np.empty([2*nx+1,2*nz+1])       #node coordinate

    for k in range(2*nz1+1):
        xg[:2*nx1,k] = np.linspace(0,10+10/(len(zg)-1)*k,2*nx1,endpoint=False)
        xg[2*nx1:2*(nx1+nx2),k] = np.linspace(10+10/(len(zg)-1)*k,area_x-(10+10/(len(zg)-1)*k),2*nx2,endpoint=False)
        xg[2*(nx1+nx2):,k] = np.linspace(area_x-(10+10/(len(zg)-1)*k),area_x,2*nx1+1,endpoint=True)      #傾斜部分

    for k in range(2*nz1+1,2*nz+1):
        xg[:,k] = np.copy(xg[:,2*nz1])



elif modelid == 2:
    area_x = 400
    area_z = 70
    dof = 2

    area_z1 = 13
    area_z2 = 40
    area_z3 = 70

    nza = 4
    nzb = 1
    nzc = 1
    nzd = 4
    nze = 1
    nzf = 1
    nzg = 1     #12-13
    nzh = 9
    nzi = 2
    nzj = 2
    nzk = 2

    nz = nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi+nzj+nzk

    nx = 200
    nxb0 = 50
    nxb1 = 55
    nxb2 = 75
    nxb3 = 20
    nxc0 = 50
    nxc1 = 55
    nxc2 = 75
    nxc3 = 20
    nxd0 = 105
    nxd1 = 75
    nxd2 = 20
    nxe0 = 105
    nxe1 = 95
    nxe2 = 20
    nxf0 = 105
    nxf1 = 95
    nxg0 = 105
    nxg1 = 95
    nxh0 = 105
    nxh1 = 95
    nxi0 = 105
    nxi1 = 95
    nxj0 = 105
    nxj1 = 95
    nxk0 = 105
    nxk1 = 95


    _zg = np.append(np.linspace(0,13,2*(nza+nzb+nzc+nzd+nze+nzf+nzg),endpoint=False),np.linspace(area_z1,area_z2,2*(nzh),endpoint=False))
    zg = np.append(_zg,np.linspace(area_z2,area_z3,2*(nzi+nzj+nzk)+1,endpoint=True))
    xg = np.empty([2*nx+1,2*nz+1])       #node coordinate

    for k in range(2*nza+1):        #0-4
        xg[:,k] = np.linspace(0,area_x,2*nx+1,endpoint=True)
    for k in range(2*nza+1,2*(nza+nzb)+1):      #4-5
        xg[:2*nxb0,k] = np.linspace(0,100-50/(2*nzb)*(k-2*nza),2*nxb0,endpoint=False)
        xg[2*nxb0:2*(nxb0+nxb1),k] = np.linspace(100-50/(2*nzb)*(k-2*nza),210,2*nxb1,endpoint=False)
        xg[2*(nxb0+nxb1):2*(nxb0+nxb1+nxb2),k] = np.linspace(210,360-10/(2*nzb)*(k-2*nza),2*nxb2,endpoint=False)
        xg[2*(nxb0+nxb1+nxb2):,k] = np.linspace(360-10/(2*nzb)*(k-2*nza),area_x,2*nxb3+1,endpoint=True)
    for k in range(2*(nza+nzb)+1,2*(nza+nzb+nzc)+1):        #5-6
        xg[:2*nxc0,k] = np.linspace(0,50+50/(2*nzc)*(k-2*(nza+nzb)),2*nxc0,endpoint=False)
        xg[2*nxc0:2*(nxc0+nxc1),k] = np.linspace(50+50/(2*nzc)*(k-2*(nza+nzb)),210,2*nxc1,endpoint=False)
        xg[2*(nxc0+nxc1):2*(nxc0+nxc1+nxc2),k] = np.linspace(210,350-10/(2*nzc)*(k-2*(nza+nzb)),2*nxc2,endpoint=False)
        xg[2*(nxc0+nxc1+nxc2):,k] = np.linspace(350-10/(2*nzc)*(k-2*(nza+nzb)),area_x,2*nxc3+1,endpoint=True)
    for k in range(2*(nza+nzb+nzc)+1,2*(nza+nzb+nzc+nzd)+1):        #6-10
        xg[:2*nxd0,k] = np.linspace(0,210+10/(2*nzd)*(k-2*(nza+nzb+nzc)),2*nxd0,endpoint=False)
        xg[2*nxd0:2*(nxd0+nxd1),k] = np.linspace(210+10/(2*nzd)*(k-2*(nza+nzb+nzc)),340-10/(2*nzd)*(k-2*(nza+nzb+nzc)),2*nxd1,endpoint=False)
        xg[2*(nxd0+nxd1):,k] = np.linspace(340-10/(2*nzd)*(k-2*(nza+nzb+nzc)),area_x,2*nxd2+1,endpoint=True)
    for k in range(2*(nza+nzb+nzc+nzd)+1,2*(nza+nzb+nzc+nzd+nze)+1):        #10-11
        xg[:2*nxe0,k] = np.linspace(0,220+120/(2*nze)*(k-2*(nza+nzb+nzc+nzd)),2*nxe0,endpoint=False)
        xg[2*nxe0:,k] = np.linspace(220+120/(2*nze)*(k-2*(nza+nzb+nzc+nzd)),area_x,2*nxe1+1,endpoint=True)
    for k in range(2*(nza+nzb+nzc+nzd+nze)+1,2*(nza+nzb+nzc+nzd+nze+nzf)+1):        #11-12
        xg[:2*nxf0,k] = np.linspace(0,340-100/(2*nzf)*(k-2*(nza+nzb+nzc+nzd+nze)),2*nxf0,endpoint=False)
        xg[2*nxf0:,k] = np.linspace(340-100/(2*nzf)*(k-2*(nza+nzb+nzc+nzd+nze)),area_x,2*nxf1+1,endpoint=True)
    for k in range(2*(nza+nzb+nzc+nzd+nze+nzf)+1,2*(nza+nzb+nzc+nzd+nze+nzf+nzg)+1):        #12-13
        xg[:2*nxg0,k] = np.linspace(0,240-90/(2*nzg)*(k-2*(nza+nzb+nzc+nzd+nze+nzf)),2*nxg0,endpoint=False)
        xg[2*nxg0:,k] = np.linspace(240-90/(2*nzg)*(k-2*(nza+nzb+nzc+nzd+nze+nzf)),area_x,2*nxg1+1,endpoint=True)
    for k in range(2*(nza+nzb+nzc+nzd+nze+nzf+nzg)+1,2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh)+1):        #13-40
        xg[:2*nxh0,k] = np.linspace(0,150+190/(2*nzh)*(k-2*(nza+nzb+nzc+nzd+nze+nzf+nzg)),2*nxh0,endpoint=False)
        xg[2*nxh0:,k] = np.linspace(150+190/(2*nzh)*(k-2*(nza+nzb+nzc+nzd+nze+nzf+nzg)),area_x,2*nxh1+1,endpoint=True)
    for k in range(2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh)+1,2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi)+1):        #40-50
        xg[:2*nxi0,k] = np.linspace(0,340-100/(2*nzi)*(k-2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh)),2*nxi0,endpoint=False)
        xg[2*nxi0:,k] = np.linspace(340-100/(2*nzi)*(k-2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh)),area_x,2*nxi1+1,endpoint=True)
    for k in range(2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi)+1,2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi+nzj)+1):        #50-60
        xg[:2*nxj0,k] = np.linspace(0,240-140/(2*nzj)*(k-2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi)),2*nxj0,endpoint=False)
        xg[2*nxj0:,k] = np.linspace(240-140/(2*nzj)*(k-2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi)),area_x,2*nxj1+1,endpoint=True)
    for k in range(2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi+nzj)+1,2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi+nzj+nzk)+1):        #50-60
        xg[:2*nxk0,k] = np.linspace(0,100+110/(2*nzk)*(k-2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi+nzj)),2*nxk0,endpoint=False)
        xg[2*nxj0:,k] = np.linspace(100+110/(2*nzk)*(k-2*(nza+nzb+nzc+nzd+nze+nzf+nzg+nzh+nzi+nzj)),area_x,2*nxk1+1,endpoint=True)


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

elif modelid == 1 or modelid == 2:
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

elif modelid == 2:
    for k in range(nz):
        im = 0
        for i in range(nx):
            im = 0
            if 4 <= k < 5 and 180 <= i < 200:
                im = 2
            if 5 <= k < 6 and (0 <= i < 50 or 180 <= i < 200):
                im = 2
            if 6 <= k < 10 and (0 <= i < 105 or 180 <= i < 200):
                im = 2
            if 10 <= k < 11:
                im = 2
            if 11 <= k:
                im = 3
                if 11 <= k < 13 and 0 <= i < 105:
                    im = 2
                if 22 <= k < 26 and 105 <= i < 200:
                    im = 1
                if 26 <= k < 28:
                    im = 1

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



nnode = inode       #number of nodes
nelem = ielem       #number of elements


### Set material ###
material_lines = []
if modelid == 0 or modelid == 1:
    material_lines += ["{} {} {} {} {} \n".format(0,"vs_vp_rho",0.0,1500.0,1750.0)]
    material_lines += ["{} {} {} {} {} \n".format(1,"vs_vp_rho",200.0,1500.0,1750.0)]

elif modelid ==2:
    material_lines += ["{} {} {} {} {} \n".format(0,"vs_vp_rho",5.0,1500.0,1750.0)]     #Fs
    material_lines += ["{} {} {} {} {} \n".format(1,"vs_vp_rho",295.0,1500.0,1750.0)]   #Kys
    material_lines += ["{} {} {} {} {} \n".format(2,"vs_vp_rho",140.0,1500.0,1750.0)]     #As
    material_lines += ["{} {} {} {} {} \n".format(3,"vs_vp_rho",135.0,1500.0,1750.0)]   #Ac,N

nmaterial = len(material_lines)


### Set output ###
output_node_lines = []
for i in range(0,nnode):
    output_node_lines += ["{} \n".format(i)]        #define output nodes

output_element_lines = []
for i in range(0,nelem-nx-len(zg)):        #define output elements
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


with open("var.in","w") as f:       #save var, depend on target area
    f.write("{} {}\n".format(modelid,"modelid"))
    f.write("{} {} {} {}\n".format(area_x,area_z,"area_x","area_z"))
    f.write("{} {} {} {}\n".format(nx,nz,"nx","nz"))
    if modelid == 1:
        f.write("{} {} {} {}\n".format(nx1,nx2,"nx1","nx2"))
        f.write("{} {} {} {}\n".format(nz1,nz2,"nz1","nz2"))
