###---USE ONLY IN "output/<YYYYMMDD-time>"---###

import numpy as np
import pandas as pd
import os
# if os.path.dirname(__file__):       #currentdirectoryをfile位置にセット
#     os.chdir(os.path.dirname(__file__))
import input.mk_mesh

##---var setup---##
modeid = 1
xnodes = input.mk_mesh.xg.shape[0]      #x-nodeがz方向に変化しない時
znodes = input.mk_mesh.zg.shape[0]
node = input.mk_mesh.node

nnodes = input.mk_mesh.output_nnode
nelem = input.mk_mesh.output_nelem




##---set nodepoints---##
points = [ "{} {} {} \n".format("POINTS",nnodes,"float") ]      #全ノード数

nodeset = []      #set node

if modeid == 0:
    for k in range(znodes):
        for i in range(xnodes):
            nodeset += [ "{} {} {} \n".format(input.mk_mesh.xg[i],0,input.mk_mesh.zg[k]) ]

elif modeid == 1:
    for k in range(znodes):
        for i in range(xnodes):
            nodeset += [ "{} {} {} \n".format(input.mk_mesh.xg[i,k],0,input.mk_mesh.zg[k]) ]


cells = [ "{} {} {} \n".format("CELLS",nelem,10*nelem) ]      #アイソ接点数+1





##---set element---##
cellset = []
if modeid == 0:
    for k in range(input.mk_mesh.nz):
        im = 1
        for i in range(input.mk_mesh.nx):
            im = 1
            if k <= 1 and i >= int(input.mk_mesh.nx/2):
                im = 0
            cellset += [ "{} {} {} {} {} {} {} {} {} {} \n".format("9",node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
                                                             node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
                                                             node[2*i+1,2*k+1])]        #9-node

elif modeid == 1:
    for k in range(input.mk_mesh.nz):
        im = 1
        for i in range(input.mk_mesh.nx):
            im = 1
            if k < input.mk_mesh.nz1 and input.mk_mesh.nx1 <= i < input.mk_mesh.nx1+input.mk_mesh.nx2:
                im = 0

            cellset += [ "{} {} {} {} {} {} {} {} {} {} \n".format("9",node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
                                                             node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
                                                             node[2*i+1,2*k+1])]        #9-node


celltypes = [ "{} {}\n".format("CELL_TYPES",nelem)]       #"CELL_TYPES,nelement"

_celltypesset = ["28 "]     #9-nodeVTK_BIQUADRATIC_QUAD
celltypesset = []
for i in range(nelem):
    celltypesset += _celltypesset


##---cell data---##
celldatasetlines1 = ["{} {}\n".format("CELL_DATA",input.mk_mesh.output_nelem)]
celldatasetlines2 = ["{} {} {}\n".format("SCALARS","strainxx","float")]





stxxdf = pd.read_table('output/output_strainxx.dat',header=None,sep="    ",usecols=lambda x: x not in[0],engine='python')
    #1列目を除き読み込み
i = 1000        #writeする行(時刻)

celldata = []
for k in range(0,len(stxxdf.columns)):
    celldata += ["{}\n".format(stxxdf.iloc[i,k])]


##--node data---##


##---write vtkfile---##
with open("vtk/output.vtk","w") as f:

    ###---mesh---###
    f.write("# vtk DataFile Version 3.0\n")
    f.write("output\n")       #vtk title
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
    f.writelines(points)        #"POINTS",nnode,dtype
    f.writelines(nodeset)       #define nodes
    f.writelines(cells)         #"CELLS",nelem,node_in_element+1
    f.writelines(cellset)       #define cells
    f.writelines(celltypes)     #"CELL_TYPES",nelem
    f.writelines(celltypesset)  #define celltype]
    f.write("\n")

    ###---cell data---###
    f.writelines(celldatasetlines1)        #"CELL_DATA",nelement
    f.writelines(celldatasetlines2)       #"SCALARS",dataname,dtype
    f.writelines("LOOKUP_TABLE default\n")      #data　reference
    f.writelines(celldata)
