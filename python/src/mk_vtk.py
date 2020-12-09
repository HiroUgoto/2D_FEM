import numpy as np
import pandas as pd
import os
if os.path.dirname(__file__):       #currentdirectoryをfile位置にセット
    os.chdir(os.path.dirname(__file__))
import input.mk_mesh

##---mesh setup---##
xnodes = str(len(input.mk_mesh.xg[0]))      #x-nodeがz方向に変化しない時
znodes = str(len(input.mk_mesh.zg))
nnodes = str(len(input.mk_mesh.xg[0])*len(input.mk_mesh.zg))

points = [ "{} {} {} \n".format("POINTS",input.mk_mesh.output_nnode,"float") ]      #全ノード数

nodeset = []      #set node
for k in range(len(input.mk_mesh.zg)):
    for i in range(len(input.mk_mesh.xg[0])):
        nodeset += [ "{} {} {} \n".format(input.mk_mesh.xg[k,i],0,input.mk_mesh.zg[k]) ]

cells = [ "{} {} {} \n".format("CELLS",input.mk_mesh.output_nelem,10*input.mk_mesh.output_nelem) ]      #アイソ接点数+1

node = input.mk_mesh.node
cellset = []
for k in range(input.mk_mesh.nz):       #set element　except1d3input
    im = 1
    for i in range(input.mk_mesh.nx):
        if k < 3:
            im = 0
            if i == 0 or i == input.mk_mesh.nx-1:
                im = 1
        cellset += [ "{} {} {} {} {} {} {} {} {} {} \n".format("9",node[2*i,2*k],node[2*i+2,2*k],node[2*i+2,2*k+2],node[2*i,2*k+2],
                                                         node[2*i+1,2*k],node[2*i+2,2*k+1],node[2*i+1,2*k+2],node[2*i,2*k+1],
                                                         node[2*i+1,2*k+1])]        #9-node
celltypes = [ "{} {}\n".format("CELL_TYPES",input.mk_mesh.output_nelem)]       #"CELL_TYPES,nelement"

_celltypesset = ["28 "]     #9-nodeVTK_BIQUADRATIC_QUAD
celltypesset = []
for i in range(input.mk_mesh.output_nelem):
    celltypesset += _celltypesset

##---cell data---##
celldatasetlines1 = ["{} {}\n".format("CELL_DATA",input.mk_mesh.output_nelem)]
celldatasetlines2 = ["{} {} {}\n".format("SCALARS","strainxx","default")]





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
    f.write("output\n")       #data title
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
    f.writelines(points)
    f.writelines(nodeset)
    f.writelines(cells)
    f.writelines(cellset)
    f.writelines(celltypes)
    f.writelines(celltypesset)
    

    ###---cell data---###
    f.writelines(celldatasetlines1)        #"CELL_DATA",nelement
    f.writelines(celldatasetlines2)       #"SCALARS",dataname,dtype
    f.writelines("LOOKUP_TABLE default\n")      #data　reference
    f.writelines(celldata)
