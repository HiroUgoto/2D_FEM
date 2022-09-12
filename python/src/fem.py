import numpy as np
import copy
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
        self.connected_elements = []

        self.enrich_elements = []
        self.normal_elements = []


    # ======================================================================= #
    def set_init(self):
        self._set_mesh()
        self._set_initial_condition()
        self._set_set()
        self._set_initial_matrix()

    # ---------------------------------------
    def _set_mesh(self):
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

            if element.material_id < 0:
                element.set_material(None)
            else:
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
            if "connect" in element.style:
                self.connected_elements += [element]

            if "X" in element.style:
                self.enrich_elements +=[element]
            else:
                self.normal_elements +=[element]

    # ---------------------------------------
    def _set_initial_condition(self):
        for node in self.nodes:
            node.set_initial_condition()

            if np.count_nonzero(node.freedom) == self.dof:
                self.free_nodes += [node]
            else:
                self.fixed_nodes += [node]

        for element in self.elements:
            element.set_pointer_list()

    # ---------------------------------------
    def _set_set(self):
        self.node_set = set(self.nodes)
        self.free_node_set = set(self.free_nodes)
        self.fixed_node_set = set(self.fixed_nodes)

        self.element_set = set(self.elements)
        self.input_element_set = set(self.input_elements)
        self.connected_element_set = set(self.connected_elements)

    # ---------------------------------------
    def _set_initial_matrix(self):
        for element in self.element_set:
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

            # if "X" in element.style:
            #     element.enrich_node_mass = np.zeros(self.dof, dtype=np.float64)
            #     element.enrich_node_c = np.zeros(self.dof, dtype=np.float64)
            #     element.enrich_node_k = np.zeros(self.dof, dtype=np.float64)
            #
            #     for i in range(self.dof):
            #         element.enrich_node_mass[i] = element.M_diag[id]
            #         element.enrich_node_c[i] = element.C_diag[id]
            #         element.enrich_node_k[i] = element.K_diag[id]
            #         id += 1

    # ======================================================================= #
    def set_output(self,outputs):
        output_node_list,output_element_list = outputs

        self.output_nnode = len(output_node_list)
        self.output_nodes = [None] * self.output_nnode
        for id,inode in enumerate(output_node_list):
            self.output_nodes[id] = self.nodes[inode]

        self.output_nelem = len(output_element_list)
        self.output_elements = [None] * self.output_nelem
        for id,ielem in enumerate(output_element_list):
            self.output_elements[id] = self.elements[ielem]

        self.output_node_set = set(self.output_nodes)
        self.output_element_set = set(self.output_elements)

    # ======================================================================= #
    def self_gravity(self):
        ### Initial condition ###
        H = 0.0
        for node in self.nodes:
            if node.xyz[1] > H:
                H = node.xyz[1]

        g,vp = 9.8,1500.0
        for node in self.node_set:
            node.u[0] = 0.0
            node.u[1] = g/(2*vp**2) * (H**2 - node.xyz[1]**2)
            node.um = np.copy(node.u)

        self._self_gravity_cg(full=False)
        self._self_gravity_cg(full=True)

        for node in self.node_set:
            node.u0 = np.copy(node.u)
            node.um = np.copy(node.u)

    # ---------------------------------------
    def _self_gravity_cg(self,full=True):
        if full:    # Calculate both components
            id = 0
        else:       # Calculate only vertical component
            id = 1

        ### CG Method ###
        for node in self.node_set:
            node.force = np.zeros(node.dof,dtype=np.float64)
        for element in self.element_set:
            element.mk_ku()
        for node in self.node_set:
            for i in range(id,node.dof):
                if node.freedom[i] == 0:
                    node._ur[i] = 0.0
                else:
                    node._ur[i] = node.static_force[i] - node.force[i]
        for element in self.connected_element_set:
            u = np.zeros(element.nodes[0].dof,dtype=np.float64)
            for node in element.node_set:
                u += node._ur
            for node in element.node_set:
                node._ur = u/element.nnode
        for element in self.input_elements:
            for node in element.node_set:
                node._ur = np.zeros(element.nodes[0].dof,dtype=np.float64)
        for node in self.node_set:
            node._up = np.copy(node._ur)
        for element in self.element_set:
            element._up = ()
            for node in element.nodes:
                element._up += (node._up.view(),)

        for it in range(10*self.nnode):
            ## y = Ap
            for node in self.node_set:
                node.force = np.zeros(node.dof,dtype=np.float64)
            for element in self.element_set:
                element.mk_ku_u(element._up)
            for node in self.node_set:
                node._uy = node.force

            ## correction boundary condition
            for node in self.node_set:
                for i in range(id,node.dof):
                    if node.freedom[i] == 0:
                        node._uy[i] = 0.0
            for element in self.connected_element_set:
                u = np.zeros(element.nodes[0].dof,dtype=np.float64)
                for node in element.node_set:
                    u += node._uy
                for node in element.node_set:
                    node._uy = u/element.nnode
            for element in self.input_elements:
                for node in element.node_set:
                    node._uy = np.zeros(element.nodes[0].dof,dtype=np.float64)

            ## alpha = rr/py
            rr,py = 0.0,0.0
            for node in self.node_set:
                rr += node._ur @ node._ur
                py += node._up @ node._uy
            alpha = rr/py

            ## x = x + alpha*p
            rr1 = 0.0
            for node in self.node_set:
                for i in range(id,node.dof):
                    if node.freedom[i] == 0:
                        pass
                    else:
                        node.u[i] += alpha*node._up[i]
                        node._ur[i] -= alpha*node._uy[i]
                rr1 += node._ur @ node._ur

            if rr1 < 1.e-10:
                break

            ## p = r + beta*p
            beta = rr1/rr
            for node in self.node_set:
                for i in range(id,node.dof):
                    if node.freedom[i] == 0:
                        pass
                    else:
                        node._up[i] = node._ur[i] + beta*node._up[i]

            if it%100 == 0:
                print(" (self gravity process .. )",it,self.nodes[0].u[1],rr1)


    # ======================================================================= #
    def update_init(self,dt):
        for node in self.node_set:
            node.inv_mc = 1.0 / (node.mass[:] + 0.5*dt*node.c[:])
            node.mass_inv_mc = node.mass[:]*node.inv_mc[:]
            node.c_inv_mc = node.c[:]*node.inv_mc[:]*0.5*dt
            node.dtdt_inv_mc = dt*dt*node.inv_mc[:]

        self.dt = dt
        self.inv_dt2 = 1./(2.*dt)
        self.inv_dtdt = 1./(dt*dt)

    # ======================================================================= #
    def update_matrix(self):
        for node in self.node_set:
            self._update_matrix_node_init(node)
        for element in self.element_set:
            self._update_matrix_set_elements(element)

    # ---------------------------------------
    def _update_matrix_node_init(self,node):
        node.mass = np.zeros(self.dof,dtype=np.float64)
        node.c    = np.zeros(self.dof,dtype=np.float64)
        node.dynamic_force = np.zeros(self.dof,dtype=np.float64)

    def _update_matrix_set_elements(self,element):
        element.set_xn()
        element.mk_local_update()

        id = 0
        for node in element.nodes:
            for i in range(self.dof):
                node.mass[i] += element.M_diag[id]
                node.c[i] += element.C_diag[id]
                node.dynamic_force[i] += element.force[id]
                id += 1

        if "X" in element.style and element.rupture:
            element.enrich_node_mass = np.zeros(self.dof*element.ncrack, dtype=np.float64)
            element.enrich_node_c = np.zeros(self.dof*element.ncrack, dtype=np.float64)
            element.enrich_node_k = np.zeros(self.dof*element.ncrack, dtype=np.float64)

            for i in range(self.dof*element.ncrack):
                element.enrich_node_mass[i] = element.M_diag[id]
                element.enrich_node_c[i] = element.C_diag[id]
                element.enrich_node_k[i] = element.K_diag[id]
                id += 1

            element.enrich_inv_mc = 1.0 / (element.enrich_node_mass[:] + 0.5*self.dt*element.enrich_node_c[:])
            element.enrich_mass_inv_mc = element.enrich_node_mass[:]*element.enrich_inv_mc[:]
            element.enrich_c_inv_mc = element.enrich_node_c[:]*element.enrich_inv_mc[:]*0.5*self.dt
            element.enrich_dtdt_inv_mc = self.dt*self.dt*element.enrich_inv_mc[:]

    def _update_matrix_set_nodes(self,node):
        node.inv_mc = 1.0 / (node.mass[:] + 0.5*self.dt*node.c[:])
        node.mass_inv_mc = node.mass[:]*node.inv_mc[:]
        node.c_inv_mc = node.c[:]*node.inv_mc[:]*0.5*self.dt
        node.dtdt_inv_mc = self.dt*self.dt*node.inv_mc[:]

    # ======================================================================= #
    def update_time(self,acc0,vel0=None,input_wave=False,FD=False):
        if FD:
            self.update_matrix()
        else:
            for node in self.node_set:
                node.dynamic_force = np.zeros(self.dof,dtype=np.float64)

        for node in self.node_set:
            self._update_time_node_init(node)

        if input_wave:
            for element in self.input_element_set:
                self._update_time_input_wave(element,vel0)
        else:
            for element in self.element_set:
                self._update_bodyforce(element,acc0)

        if FD:
            for element in self.element_set:
                element.mk_B_stress()
                element.mk_cv()
        else:
            for element in self.element_set:
                element.mk_ku_cv()

        for node in self.free_node_set:
            self._update_time_set_free_nodes(node)
        for node in self.fixed_node_set:
            self._update_time_set_fixed_nodes(node)

        for element in self.connected_element_set:
            self._update_time_set_connected_elements_(element)

        for element in self.output_element_set:
            element.calc_stress()

    # ======================================================================= #
    def update_time_disp(self,forced_disp0,forced_nodes):
        for element in self.enrich_elements:
            element.check_enrich_rupture()

        self.update_matrix()

        for node in self.node_set:
            node.dynamic_force = np.zeros(self.dof,dtype=np.float64)
            self._update_time_node_init(node)

        for element in self.normal_elements:
            element.mk_ku_cv()
        for element in self.enrich_elements:
            element.mk_ku_cv_enrich()

        for node in self.free_node_set:
            self._update_time_set_free_nodes(node)
        for node in self.fixed_node_set:
            self._update_time_set_fixed_nodes(node)

        for element in self.connected_element_set:
            self._update_time_set_connected_elements_(element)

        for element in self.enrich_elements:
            self._update_time_set_enrich_nodes(element)

        for id in forced_nodes:
            for i in range(self.dof):
                if self.nodes[id].freedom[i] == 0:
                    self.nodes[id].u[i] = forced_disp0[i]

        for element in self.output_element_set:
            element.calc_stress()


    # ---------------------------------------
    def _update_time_node_init(self,node):
        node.force = -node.dynamic_force.copy()

    def _update_time_input_wave(self,element,vel0):
        cv = element.C @ np.tile(vel0,element.nnode)
        for i in range(element.nnode):
            i0 = self.dof*i
            element.nodes[i].force[:] -= 2*cv[i0:i0+self.dof]

    def _update_bodyforce(self,element,acc0):
        element.mk_bodyforce(acc0)
        for i in range(element.nnode):
            i0 = self.dof*i
            element.nodes[i].force[:] -= element.force[i0:i0+self.dof]

    def _update_time_set_free_nodes(self,node):
        u = np.copy(node.u)
        node.u[:] = node.mass_inv_mc*(2.*u-node.um) + node.c_inv_mc*node.um - node.dtdt_inv_mc*node.force
        node.v[:] = (node.u - node.um) * self.inv_dt2
        node.a[:] = (node.u - 2.*u + node.um) * self.inv_dtdt
        node.um = u

    def _update_time_set_fixed_nodes(self,node):
        u = np.copy(node.u)
        for i in range(self.dof):
            if node.freedom[i] == 0:
                node.u[i] = 0.0
                node.v[i] = 0.0
                node.a[i] = 0.0
            else:
                node.u[i] = node.mass_inv_mc[i]*(2.*u[i]-node.um[i]) + node.c_inv_mc[i]*node.um[i] - node.dtdt_inv_mc[i]*node.force[i]
                node.v[i] = (node.u[i] - node.um[i]) * self.inv_dt2
                node.a[i] = (node.u[i] - 2.*u[i] + node.um[i]) * self.inv_dtdt
        node.um = u

    def _update_time_set_connected_elements_(self,element):
        u = np.zeros_like(element.nodes[0].u)
        a = np.zeros_like(element.nodes[0].a)
        for node in element.node_set:
            u[:] += node.u[:]
            a[:] += node.a[:]
        for node in element.node_set:
            node.u[:] = u[:]/element.nnode
            node.a[:] = a[:]/element.nnode

    # ------------------------------------------------
    def _update_time_set_enrich_nodes(self,element):
        if element.rupture:

            id = 0
            for ic in range(element.ncrack):
                for i in range(self.dof):
                    du = np.copy(element.du_list[ic][i])

                    element.du_list[ic][i] = element.enrich_mass_inv_mc[id]*(2.*du-element.dum_list[ic][i]) + element.enrich_c_inv_mc[id]*element.dum_list[ic][i] - element.enrich_dtdt_inv_mc[id]*element.enrich_force[id]
                    element.dv_list[ic][i] = (element.du_list[ic][i] - element.dum_list[ic][i]) * self.inv_dt2
                    id += 1

                    element.dum_list[ic][i] = np.copy(du)


    # ======================================================================= #
    def print_all(self):
        for node in self.nodes:
            node.print()
        for element in self.elements:
            element.print()
