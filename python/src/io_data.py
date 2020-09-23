import numpy as np
import fem
import node,element


def input_mesh(mesh_file):
    with open(mesh_file) as f:
        lines = f.readlines()

        nnode,nelem,dof = [int(s) for s in lines[0].split()]

        irec = 1
        nodes = []
        for inode in range(nnode):
            items = lines[inode+irec].split()

            id = int(items[0])
            xyz = [float(s) for s in items[1:3]]
            freedom = [int(s) for s in items[3:]]

            nodes += [node.Node(id,xyz,freedom)]

        irec += nnode
        elements = []
        for ielem in range(nelem):
            items = lines[ielem+irec].split()

            id = int(items[0])
            style = int(items[1])
            param = [float(s) for s in items[2:5]]
            inode = [int(s) for s in items[5:]]

            elements += [element.Element(id,style,param,inode)]

        return fem.Fem(dof,nodes,elements)
