import numpy as np
import matplotlib.pyplot as plt
import os


modelid = 0

if modelid == 0:       #容器モデル
    area_x = 5.0
    area_z = 5.0
    nx =  5
    nz =  5
    dof = 2

    xg = np.linspace(0,area_x,2*nx+1,endpoint=True)
    zg = np.linspace(0,area_z,2*nz+1,endpoint=True)

elif modelid == 1:      #凹型モデル
    area_x = 12.0
    area_z = 6.0
    nx1,nx2 = 3, 6
    nz1,nz2,nz3 = 2, 3, 1
    nx =  nx1 + nx2 + nx1
    nz =  nz1 + nz2 + nz3
    dof = 2

    xg = np.linspace(0,area_x,2*nx+1,endpoint=True)
    zg = np.linspace(0,area_z,2*nz+1,endpoint=True)


### Set node ###
node = np.empty([len(xg),len(zg)],dtype=np.int32)       #node_idを振った配列(転置)
node_lines = []
plot_node = []

inode = 0
if modelid == 0:
    for k in range(len(zg)):
        for i in range(len(xg)):
            dofx,dofz = 1,1
            if k == len(zg)-1:
                dofz = 0
            if i == 0:
                dofx = 0
            if i == len(xg)-1:
                dofx = 0
            
            plot_node.append([inode,xg[i],zg[k]])

            node[i,k] = inode
            node_lines += [ "{} {} {} {} {}\n".format(inode,xg[i],zg[k],dofx,dofz)]
            inode += 1

elif modelid == 1:
    node = np.empty([len(xg),len(zg)],dtype=np.int32)       #node_idを振った配列(転置)
    node_lines = []
    for k in range(len(zg)):
        for i in range(len(xg)):
            dofx,dofz = 1,1
            if k == len(zg)-1:
                dofz = 0

            plot_node.append([inode,xg[i],zg[k]]) 

            node[i,k] = inode
            node_lines += [ "{} {} {} {} {}\n".format(inode,xg[i],zg[k],dofx,dofz) ]          
            inode += 1


### Set element ###
element_lines = []
plot_line = []


ielem = 0
if modelid == 0:
    for k in range(nz):
        for i in range(nx):
            im = 0

            style = "2d9solid"

            param_line = "{} {} {} ".format(ielem,style,im)
            style_line = "{} {} {} {} {} {} {} {} {}".format(node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
                                                            node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
                                                            node[2*i+1,2*k+1])

            element_lines += [param_line + style_line + "\n"]

            plot_line.append([im,plot_node[node[2*i,2*k]][1],plot_node[node[2*i,2*k]][2],
                        plot_node[node[2*i+2,2*k]][1],plot_node[node[2*i+2,2*k]][2],
                        plot_node[node[2*i+2,2*k+2]][1],plot_node[node[2*i+2,2*k+2]][2],
                        plot_node[node[2*i,2*k+2]][1],plot_node[node[2*i,2*k+2]][2]])

            ielem += 1

elif modelid == 1:
    for k in range(nz):
        for i in range(nx):
            im = 2
            if k < 1:
                im = 1
            elif 1 <= k < 2 or (2 <= k < 5 and nx1 <= i < nx1+nx2):
                im = 0

            style = "2d9solid"

            param_line = "{} {} {} ".format(ielem,style,im)
            style_line = "{} {} {} {} {} {} {} {} {}".format(node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
                                                                node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
                                                                node[2*i+1,2*k+1])

            element_lines += [param_line + style_line + "\n"]
            
            plot_line.append([im,plot_node[node[2*i,2*k]][1],plot_node[node[2*i,2*k]][2],
                                    plot_node[node[2*i+2,2*k]][1],plot_node[node[2*i+2,2*k]][2],
                                    plot_node[node[2*i+2,2*k+2]][1],plot_node[node[2*i+2,2*k+2]][2],
                                    plot_node[node[2*i,2*k+2]][1],plot_node[node[2*i,2*k+2]][2]])

            ielem += 1

    # for i in range(nx):
    #     style = "1d3input"
    #     im = 0

    #     param_line = "{} {} {} ".format(ielem,style,im)
    #     style_line = "{} {} {}".format(node[2*i,-1],node[2*i+2,-1],node[2*i+1,-1])

    #     element_lines += [param_line + style_line + "\n"]
    #     ielem += 1

    # for k in range(len(zg)):     #connected element
    #     style = "connect"
    #     im = -1

    #     param_line = "{} {} {} ".format(ielem,style,im)
    #     style_line = "{} {}".format(node[0,k],node[2*nx,k])

    #     element_lines += [param_line + style_line + "\n"]
    #     ielem += 1


nnode = inode       #number of nodes
nelem = ielem       #number of elements


### Set material ###
material_lines = []
if modelid == 0:
    material_lines += ["{} {} {} {} {} \n".format(0,"vs_vp_rho",10.0,1500.0,1000.0)]

elif modelid == 1:
    material_lines += ["{} {} {} {} {} \n".format(0,"vs_vp_rho",0.0,1500.0,1836.7)]     #液状化層
    material_lines += ["{} {} {} {} {} \n".format(1,"vs_vp_rho",140.0,279.8,1734.7)]   #盛土層
    material_lines += ["{} {} {} {} {} \n".format(2,"vs_vp_rho",140.0,1500.0,1836.7)]   #非液状化層

nmaterial = len(material_lines)

### Set output ###
output_node_lines = []
for i in range(0,nnode):
    output_node_lines += ["{}\n".format(i)]        #define output nodes

output_element_lines = []
for i in range(0,nelem-nx-len(zg)):        #define output elements
    output_element_lines += ["{} \n".format(i)]

output_nnode = len(output_node_lines)
output_nelem = len(output_element_lines)


### write file ###
with open("mesh.in","w",newline="\n") as f:
    f.write("{} {} {} {}\n".format(nnode,nelem,nmaterial,dof))
    f.writelines(node_lines)
    f.writelines(element_lines)
    f.writelines(material_lines)

with open("output.in","w",newline="\n") as f:
    f.write("{} {}\n".format(output_nnode,output_nelem))
    f.writelines(output_node_lines)
    f.writelines(output_element_lines)

# with open("var.in","w",newline="\n") as f:       #for mk_vtks
#     f.write("{} {}\n".format(modelid,"modelid"))
#     f.write("{} {} {} {}\n".format(area_x,area_z,"area_x","area_z"))
#     f.write("{} {} {} {}\n".format(nx,nz,"nx","nz"))
#     if modelid == 1:
#         f.write("{} {} {} {}\n".format(nx1,nx2,"nx1","nx2"))
#         f.write("{} {} {} {}\n".format(nz1,nz2,"nz1","nz2"))


### plot model ###
def plot_mesh(amp=1.0):
    pc = ["lightblue","gray","yellow","green","pink"]

    fig,ax = plt.subplots(figsize=(6,4))

    # x = [node.xyz[0] for node in fem.nodes]
    # z = [node.xyz[1] for node in fem.nodes]

    ax.set_xlim([min(xg)-0.1*area_x,max(xg)+0.1*area_x])
    ax.set_ylim([max(zg)+0.1*area_z,min(zg)-0.25*area_z])
    ax.set_aspect('equal')

    for i in range(len(plot_line)):
        ic = plot_line[i][0] % len(pc)

        f0 = (plot_line[i][1], plot_line[i][2])
        f1 = (plot_line[i][3], plot_line[i][4])
        f2 = (plot_line[i][5], plot_line[i][6])
        f3 = (plot_line[i][7], plot_line[i][8])

        fpoly = plt.Polygon((f0,f1,f2,f3),ec="k",fc=pc[ic])
        ax.add_patch(fpoly)

    rc = 0.01*area_z
    for i in range(len(plot_node)):
        p = plt.Circle((plot_node[i][1],plot_node[i][2]),rc,color="k")
        ax.add_patch(p)

    plt.show()

plot_mesh()