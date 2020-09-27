import numpy as np

class Fem():
    def __init__(self,dof,nodes,elements):
        self.nnode = len(nodes)
        self.nelem = len(elements)
        self.dof = dof

        self.nodes = nodes
        self.elements = elements

    def set_init(self):
        self.set_mesh()
        self.set_initial_condition()
        self.set_initial_matrix()

    # ------------------------------------------------
    def set_mesh(self):
        for element in self.elements:
            nodes = []
            for inode in element.inode:
                if self.nodes[inode].id == inode:
                    nodes += [self.nodes[inode]]
                else:
                    for n in self.nodes:
                        if n.id == inode:
                            nodes += [n]
                            break

            element.set_nodes(nodes)

    def set_initial_condition(self):
        for node in self.nodes:
            node.set_initial_condition()

    def set_initial_matrix(self):
        for element in self.elements:
            element.mk_local_matrix(self.dof)

            id = 0
            for node in element.nodes:
                for i in range(self.dof):
                    node.mass[i] += element.M_diag[id]
                    node.c[i] += element.C_diag[id]
                    id += 1

    # ------------------------------------------------
    def update_init(self,dt):
        for node in self.nodes:
            node.inv_mc[:] = 1.0 / (node.mass[:] + 0.5*dt*node.c[:])


    def update_time(self,dt,acc0,vel0=0.0,input_wave=False):
        for node in self.nodes:
            node.force = np.zeros(self.dof,dtype=np.float64)

        if input_wave:
            for element in self.elements:
                v = np.zeros([element.dof*element.nnode])
                for i in range(element.nnode):
                    i0 = self.dof*i
                    v[i0] = vel0

                cv = np.dot(element.C+np.diag(element.C_diag),v)
                for i in range(element.nnode):
                    i0 = self.dof*i
                    element.nodes[i].force[:] -= 2*cv[i0:i0+self.dof]

        else:
            for node in self.nodes:
                node.force += np.dot(np.diag(node.mass),np.array([acc0,0.0]))

        for element in self.elements:
            element.mk_ku(self.dof)
            element.mk_cv(self.dof)

        for node in self.nodes:
            for i in range(self.dof):
                if node.freedom[i] == 0:
                    node.u[i]  = 0.0
                    node.um[i] = 0.0
                    node.v[i]  = 0.0

                else:
                    um = np.copy(node.um[i])
                    u = np.copy(node.u[i])

                    node.u[i] = (node.mass[i]*(2*u-um) + 0.5*dt*node.c[i]*um - dt*dt*node.force[i]) * node.inv_mc[i]

                    node.v[i] = (node.u[i] - um) / (2*dt)
                    node.um[i] = np.copy(u)



    # ------------------------------------------------
    def print_all(self):
        for node in self.nodes:
            node.print()
        for element in self.elements:
            element.print()
