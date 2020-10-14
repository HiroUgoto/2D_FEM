import numpy as np
import fem
import node,element,material

# ------------------------------------------------------------------- #
def input_mesh(mesh_file):
    with open(mesh_file) as f:
        lines = f.readlines()

        nnode,nelem,nmaterial,dof = [int(s) for s in lines[0].split()]

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
            style = items[1]
            material_id = int(items[2])
            inode = [int(s) for s in items[3:]]

            elements += [element.Element(id,style,material_id,inode)]

        irec += nelem
        materials = []
        for imaterial in range(nmaterial):
            items = lines[imaterial+irec].split()

            id = int(items[0])
            style = items[1]
            param = [float(s) for s in items[2:]]

            materials += [material.Material(id,style,param)]

        return fem.Fem(dof,nodes,elements,materials)

# ------------------------------------------------------------------- #
def input_outputs(output_file):
    with open(output_file) as f:
        lines = f.readlines()

        nnode,nelem = [int(s) for s in lines[0].split()]

        irec = 1
        nodes = []
        for inode in range(nnode):
            items = lines[inode+irec].split()

            inode = int(items[0])
            nodes += [inode]

        irec += nnode
        elements = []
        for ielem in range(nelem):
            items = lines[ielem+irec].split()

            ielem = int(items[0])
            elements += [ielem]

        return nodes, elements
