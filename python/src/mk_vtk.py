###---USE ONLY IN "output/<YYYYMMDD-time>"---###

import numpy as np
import pandas as pd
import os
# if os.path.dirname(__file__):       #currentdirectoryをfile位置にセット
#     os.chdir(os.path.dirname(__file__))

##--read files--##
with open("var.in") as f:
    lines = f.readlines()
    modelid = lines[0][0]
    area_x = lines[1][0]
    area_z = lines[1][1]
    nx = lines[2][0]
    nz = lines[2][1]
    if modelid == 1:
        nx1 = lines[3][0]
        nx2 = lines[3][1]
        nz1 = lines[4][0]
        nz2 = lines[4][0]


with open("mesh.in") as f:
    lines = f.readlines()
    nnode,nelem,nmaterial,dof = [int(s) for s in lines[0].split()]  #cautin nelem

    irec = 1
    nodes = []
    for inode in range(nnode):
        items = lines[inode+irec].split()
        nodes += [items]        #[id,x,z,dofx,dofz]

    irec += nnode
    elements = []
    for ielem in range(nelem):
        items = lines[ielem+irec].split()
        elements += [items]         #[id,style,im,nodeid


    irec += nelem
    materials = []
    for imaterial in range(nmaterial):
        items = lines[imaterial+irec].split()

        id = int(items[0])
        style = items[1]
        param = [float(s) for s in items[2:11]]




##---var setup---##





##---set nodepoints---##
points = [ "{} {} {} \n".format("POINTS",nnode,"float") ]      #全ノード数

nodeset = []      #set node
for i in range(nnode):
    nodeset += [ "{} {} {}\n".format(nodes[i][1],0,nodes[i][2])]

##---set element---##
cellset = []
for i in range(nelem):
    if elements[i][1] == "2d9solid":
        cellset += [ "{} {}\n".format("9"," ".join(elements[i][3:]))]

nselem = len(cellset)
cells = [ "{} {} {} \n".format("CELLS",nselem,10*nselem) ]      #アイソ接点数+1

celltypes = [ "{} {}\n".format("CELL_TYPES",nselem)]       #"CELL_TYPES,nelement"

_celltypesset = ["28 "]     #9-nodeVTK_BIQUADRATIC_QUAD
celltypesset = []
for i in range(nselem):
    celltypesset += _celltypesset


##---cell data---##
celldatasetlines1 = ["{} {}\n".format("CELL_DATA",nselem)]
celldatasetlines2 = ["{} {} {}\n".format("SCALARS","strainxx","float")]





stxxdf = pd.read_table('strainxx.dat',header=None,sep="    ",usecols=lambda x: x not in[0],engine='python')
    #1列目を除き読み込み
i = 100        #writeする行(時刻)

celldata = []
for k in range(0,len(stxxdf.columns)):
    celldata += ["{}\n".format(stxxdf.iloc[i,k])]


##--node data---##


##---write vtkfile---##
with open("output.vtk","w") as f:

    ###---mesh---###
    f.write("# vtk DataFile Version 3.0\n")
    f.write("output\n")       #vtk title
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
    f.writelines(points)        #"POINTS",nnode,dtype
    f.writelines(nodeset)       #define nodes
    f.writelines(cells)         #"CELLS",nselem,node_in_element+1
    f.writelines(cellset)       #define cells
    f.writelines(celltypes)     #"CELL_TYPES",nselem
    f.writelines(celltypesset)  #define celltype]
    f.write("\n")

    ###---cell data---###
    f.writelines(celldatasetlines1)        #"CELL_DATA",nselem
    f.writelines(celldatasetlines2)       #"SCALARS",dataname,dtype
    f.writelines("LOOKUP_TABLE default\n")      #data　reference
    f.writelines(celldata)
