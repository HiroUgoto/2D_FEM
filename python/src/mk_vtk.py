import numpy as np
import input.mk_mesh


xnodes = str(len(input.mk_mesh.xg))
znodes = str(len(input.mk_mesh.zg))
dimensions = [ "{} {} {} {} \n".format("DIMENSIONS",xnodes,"1",znodes) ]        #いらんくね？
points = [ "{} {} {} \n".format("POINTS",input.mk_mesh.output_nnode,"float") ]      #全ノード数
nodeset = []      #set node
for k in range(len(input.mk_mesh.zg)):
    for i in range(len(input.mk_mesh.xg)):
        nodeset += [ "{} {} {} \n".format(input.mk_mesh.xg[i],0,input.mk_mesh.zg[k]) ]
cells = [ "{} {} {} \n".format("CELLS",input.mk_mesh.output_nelem,"10") ]      #アイソ接点数+1

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



with open("vtk\output.vtk","w") as f:
    f.write("# vtk DataFile Version 3.0\n")
    f.write("output\n")       #data title
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
#    f.writelines(dimensions)
    f.writelines(points)
    f.writelines(nodeset)
    f.writelines(cells)
    f.writelines(cellset)
