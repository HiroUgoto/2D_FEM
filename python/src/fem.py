import numpy as np
from concurrent import futures

class Fem():
    def __init__(self,dof,nodes,elements,materials):
        self.nnode = len(nodes)
        self.nelem = len(elements)
        self.dof = dof

        self.nodes = nodes
        self.elements = elements
        self.materials = materials

        self.input_elements = []
        self.free_nodes = []
        self.fixed_nodes = []

    def set_init(self):
        self.set_mesh()
        self.set_initial_condition()
        self.set_initial_matrix()

    # ------------------------------------------------
    def set_output(self,outputs):
        output_node_list,output_element_list = outputs

        self.output_nodes = []
        for inode in output_node_list:
            self.output_nodes += [self.nodes[inode]]
        self.output_nnode = len(self.output_nodes)

        self.output_elements = []
        for ielem in output_element_list:
            self.output_elements += [self.elements[ielem]]
        self.output_nelem = len(self.output_elements)

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

            if self.materials[element.material_id].id == element.material_id:
                material = self.materials[element.material_id]
            else:
                for m in self.materials:
                    if m.id == element.material_id:
                        material = m
                        break
            element.set_material(material)

            if "input" in element.style:
                self.input_elements += [element]

    # ------------------------------------------------
    def set_initial_condition(self):
        for node in self.nodes:
            node.set_initial_condition()

            if np.count_nonzero(node.freedom) == self.dof:
                self.free_nodes += [node]
            else:
                self.fixed_nodes += [node]

        for element in self.elements:
            element.set_pointer_list()

    # ------------------------------------------------
    def set_initial_matrix(self):
        for element in self.elements:
            element.set_xn()
            element.mk_local_matrix_init(self.dof)
            element.mk_local_matrix()
            element.mk_local_vector()

            id = 0
            for node in element.nodes:
                for i in range(self.dof):
                    node.mass[i] += element.M_diag[id]
                    node.c[i] += element.C_diag[id]
                    node.k[i] += element.K_diag[id]
                    node.static_force[i] += element.force[id]
                    id += 1

    # ------------------------------------------------
    # ------------------------------------------------
    def self_gravity(self):
        ### Initial condition ###
        H = 0.0
        for node in self.nodes:
            if node.xyz[1] > H:
                H = node.xyz[1]

        g,vp = 9.8,1500.0
        for node in self.nodes:
            node.u[0] = 0.0
            node.u[1] = g/(2*vp**2) * (H**2 - node.xyz[1]**2)
            node.um = np.copy(node.u)

        ### CG Method ###
        for node in self.nodes:
            node.force = np.zeros(node.dof,dtype=np.float64)
        for element in self.elements:
            element.mk_ku()
        for node in self.nodes:
            for i in range(1,node.dof):
                if node.freedom[i] == 0:
                    node._ur[i] = 0.0
                else:
                    node._ur[i] = node.static_force[i] - node.force[i]
                node._up = np.copy(node._ur)
        for element in self.elements:
            element._up = ()
            for node in element.nodes:
                element._up += (node._up.view(),)

        for it in range(10*self.nnode):
            ## y = Ap
            for node in self.nodes:
                node._uy = np.zeros(node.dof,dtype=np.float64)
            for element in self.elements:
                ku = np.dot(element.K,np.hstack(element._up))
                for i in range(element.nnode):
                    i0 = element.dof*i
                    element.nodes[i]._uy[:] += ku[i0:i0+element.dof]

            ## alpha = rr/py
            rr,py = 0.0,0.0
            for node in self.nodes:
                for i in range(1,node.dof):
                    if node.freedom[i] == 0:
                        node._uy[i] = 0.0
                rr += np.dot(node._ur,node._ur)
                py += np.dot(node._up,node._uy)
            alpha = rr/py

            ## x = x + alpha*p
            rr1 = 0.0
            for node in self.nodes:
                for i in range(1,node.dof):
                    if node.freedom[i] == 0:
                        pass
                    else:
                        node.u[i] += alpha*node._up[i]
                        node._ur[i] -= alpha*node._uy[i]
                rr1 += np.dot(node._ur,node._ur)

            if rr1 < 1.e-10:
                break

            ## p = r + beta*p
            beta = rr1/rr
            for node in self.nodes:
                for i in range(1,node.dof):
                    if node.freedom[i] == 0:
                        pass
                    else:
                        node._up[i] = node._ur[i] + beta*node._up[i]

            # if it%100 == 0:
            print(self.nodes[0].u[1],rr1)


        for node in self.nodes:
            node.u0 = np.copy(node.u)


    # ------------------------------------------------
    def update_init(self,dt):
        for node in self.nodes:
            node.inv_mc = 1.0 / (node.mass[:] + 0.5*dt*node.c[:])
            node.mass_inv_mc = node.mass[:]*node.inv_mc[:]
            node.c_inv_mc = node.c[:]*node.inv_mc[:]*0.5*dt
            node.dtdt_inv_mc = dt*dt*node.inv_mc[:]

        self.dt = dt
        self.inv_dt2 = 1./(2.*dt)

    # ------------------------------------------------
    # ------------------------------------------------
    def update_matrix(self):
        for node in self.nodes:
            self._update_matrix_node_init(node)
        for element in self.elements:
            self._update_matrix_set_elements(element)
        for node in self.nodes:
            self._update_matrix_set_nodes(node)

    # ------------------------------------------------
    def _update_matrix_node_init(self,node):
        node.mass = np.zeros(self.dof,dtype=np.float64)
        node.c    = np.zeros(self.dof,dtype=np.float64)
        node.dynamic_force = np.zeros(self.dof,dtype=np.float64)
        # node.dynamic_force = -np.copy(node.static_force)

    def _update_matrix_set_elements(self,element):
        element.set_xn()
        element.mk_local_matrix()
        element.mk_local_vector()

        id = 0
        for node in element.nodes:
            for i in range(self.dof):
                node.mass[i] += element.M_diag[id]
                node.c[i] += element.C_diag[id]
                node.dynamic_force[i] += element.force[id]
                id += 1

    def _update_matrix_set_nodes(self,node):
        node.inv_mc = 1.0 / (node.mass[:] + 0.5*self.dt*node.c[:])
        node.mass_inv_mc = node.mass[:]*node.inv_mc[:]
        node.c_inv_mc = node.c[:]*node.inv_mc[:]*0.5*self.dt
        node.dtdt_inv_mc = self.dt*self.dt*node.inv_mc[:]

    # ------------------------------------------------
    # ------------------------------------------------
    def update_time(self,acc0,vel0=None,input_wave=False,FD=False):
        if FD:
            self.update_matrix()   # Finite deformation
        else:
            for node in self.nodes:
                node.dynamic_force = np.copy(node.static_force)

        for node in self.nodes:
            self._update_time_node_init(node)

        if input_wave:
            for element in self.input_elements:
                self._update_time_input_wave(element,vel0)
        else:
            for element in self.elements:
                self._update_bodyforce(element,acc0)
            # for node in self.nodes:
            #     node.force += np.dot(np.diag(node.mass),acc0)

        if FD:
            for element in self.elements:
                element.mk_B_stress()    # Finite deformation
                element.mk_cv()
        else:
            for element in self.elements:
                element.mk_ku()        # Not finite deformation
                element.mk_cv()

        for node in self.free_nodes:
            self._update_time_set_free_nodes(node)
        for node in self.fixed_nodes:
            self._update_time_set_fixed_nodes(node)

        for element in self.output_elements:
            element.calc_stress()

    # ------------------------------------------------
    def _update_time_node_init(self,node):
        node.force = -node.dynamic_force.copy()

    def _update_time_input_wave(self,element,vel0):
        cv = np.dot(element.C,np.tile(vel0,element.nnode))
        for i in range(element.nnode):
            i0 = self.dof*i
            element.nodes[i].force[:] -= 2*cv[i0:i0+self.dof]

    def _update_bodyforce(self,element,acc0):
        element.mk_bodyforce(acc0)
        for i in range(element.nnode):
            i0 = self.dof*i
            element.nodes[i].force[:] -= element.force[i0:i0+self.dof]

    def _update_time_set_free_nodes(self,node):
        u = node.u.copy()
        node.u[:] = node.mass_inv_mc*(2.*u-node.um) + node.c_inv_mc*node.um - node.dtdt_inv_mc*node.force
        node.v[:] = (node.u - node.um) * self.inv_dt2
        node.um = u

    def _update_time_set_fixed_nodes(self,node):
        u = node.u.copy()
        for i in range(self.dof):
            if node.freedom[i] == 0:
                node.u[i]  = 0.0
            else:
                node.u[i] = node.mass_inv_mc[i]*(2.*u[i]-node.um[i]) + node.c_inv_mc[i]*node.um[i] - node.dtdt_inv_mc[i]*node.force[i]
        node.v[:] = (node.u - node.um) * self.inv_dt2
        node.um = u


    # ------------------------------------------------
    # ------------------------------------------------
    def print_all(self):
        for node in self.nodes:
            node.print()
        for element in self.elements:
            element.print()
