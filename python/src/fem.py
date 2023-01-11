import numpy as np
import concurrent.futures
import sys

class Fem():
    def __init__(self,dof,nodes,elements,materials):
        self.nnode = len(nodes)
        self.nelem = len(elements)
        self.dof = dof

        self.nodes = nodes
        self.elements = elements
        self.materials = materials

        self.free_nodes = []
        self.fixed_nodes = []

        self.input_elements = []
        self.connected_elements = []
        self.spring_elements = []

        self.e_elements = []
        self.ep_elements = []
        self.ep_eff_elements = []

    # ======================================================================= #
    def set_init(self):
        self._set_mesh()
        self._set_initial_condition()
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
                element.set_material(self.dof,None)
            else:
                if self.materials[element.material_id].id == element.material_id:
                    material = self.materials[element.material_id]
                else:
                    for m in self.materials:
                        if m.id == element.material_id:
                            material = m
                            break

                element.set_material(self.dof,material)

                if "ep_eff_" in material.style:
                    self.ep_eff_elements += [element]
                elif "ep_" in material.style:
                    self.ep_elements += [element]
                else:
                    self.e_elements += [element]

            if "input" in element.style:
                self.input_elements += [element]
            if "connect" in element.style:
                self.connected_elements += [element]
            if "spring" in element.style:
                self.spring_elements += [element]

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
    def _set_initial_matrix(self):
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

    # ======================================================================= #
    def set_output(self,outputs):
        output_node_list,output_element_list = outputs

        self.output_nnode = len(output_node_list)
        self.output_nodes = [None] * self.output_nnode
        for id,inode in enumerate(output_node_list):
            self.output_nodes[id] = self.nodes[inode]

        self.output_nelem = len(output_element_list)
        self.output_elements = [None] * self.output_nelem
        self.output_e_elements,self.output_ep_elements = [],[]
        for id,ielem in enumerate(output_element_list):
            self.output_elements[id] = self.elements[ielem]
            if "ep_" in self.elements[ielem].material.style:
                self.output_ep_elements += [self.elements[ielem]]
            else:
                self.output_e_elements += [self.elements[ielem]]

    # ======================================================================= #
    def set_rayleigh_damping(self,f0,f1,h):
        omega0,omega1 = 2*np.pi*f0, 2*np.pi*f1
        alpha = 2*h*omega0*omega1 / (omega0+omega1)
        beta = 2*h / (omega0+omega1)
        print(alpha,beta)

        for node in self.nodes:
            node.c    = np.zeros(self.dof,dtype=np.float64)

        for element in self.elements:
            element.mk_local_damping_matrix(alpha,beta)

            id = 0
            for node in element.nodes:
                for i in range(self.dof):
                    node.c[i] += element.C_diag[id]
                    id += 1


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
    def _self_gravity_cg(self,full=True,print_flag=True):
        if full:    # Calculate both components
            id = 0
        else:       # Calculate only vertical component
            id = 1

        ### CG Method ###
        for node in self.nodes:
            node.force = np.zeros(node.dof,dtype=np.float64)
        for element in self.elements:
            element.mk_ku()
        for node in self.nodes:
            for i in range(id,node.dof):
                if node.freedom[i] == 0:
                    node._ur[i] = 0.0
                else:
                    node._ur[i] = node.static_force[i] - node.force[i]
        for element in self.connected_elements:
            u = np.zeros(element.nodes[0].dof,dtype=np.float64)
            for node in element.nodes:
                u += node._ur
            for node in element.nodes:
                node._ur = u/element.nnode
        for element in self.input_elements:
            for node in element.nodes:
                node._ur = np.zeros(element.nodes[0].dof,dtype=np.float64)
        for node in self.nodes:
            node._up = np.copy(node._ur)
        for element in self.elements:
            element._up = ()
            for node in element.nodes:
                element._up += (node._up.view(),)

        for it in range(10*self.nnode):
            ## y = Ap
            for node in self.nodes:
                node.force = np.zeros(node.dof,dtype=np.float64)
            for element in self.elements:
                element.mk_ku_u(element._up)
            for node in self.nodes:
                node._uy = node.force

            ## correction boundary condition
            for node in self.nodes:
                for i in range(id,node.dof):
                    if node.freedom[i] == 0:
                        node._uy[i] = 0.0
            for element in self.connected_elements:
                u = np.zeros(element.nodes[0].dof,dtype=np.float64)
                for node in element.nodes:
                    u += node._uy
                for node in element.nodes:
                    node._uy = u/element.nnode
            for element in self.input_elements:
                for node in element.nodes:
                    node._uy = np.zeros(element.nodes[0].dof,dtype=np.float64)


            ## alpha = rr/py
            rr,py = 0.0,0.0
            for node in self.nodes:
                rr += node._ur @ node._ur
                py += node._up @ node._uy
            alpha = rr/py

            ## x = x + alpha*p
            rr1 = 0.0
            for node in self.nodes:
                for i in range(id,node.dof):
                    if node.freedom[i] == 0:
                        pass
                    else:
                        node.u[i] += alpha*node._up[i]
                        node._ur[i] -= alpha*node._uy[i]
                rr1 += node._ur @ node._ur

            if rr1 < 1.e-10:
                if print_flag:
                    print(" (self gravity process .. )",it,self.nodes[0].u[1],rr1)
                break

            ## p = r + beta*p
            beta = rr1/rr
            for node in self.nodes:
                for i in range(id,node.dof):
                    if node.freedom[i] == 0:
                        pass
                    else:
                        node._up[i] = node._ur[i] + beta*node._up[i]

            if it%100 == 0 and print_flag:
                print(" (self gravity process .. )",it,self.nodes[0].u[1],rr1)

    # ======================================================================= #
    def set_ep_initial_state(self):
        self._self_gravity_cg(full=False,print_flag=False)

        for element in self.ep_elements:
            element.calc_stress()
            # element.ep.initial_state_isotropic(element.stress)
            element.ep.initial_state(element.stress)
            element.material.rmu,element.material.rlambda = element.ep.elastic_modulus()
            element.clear_strain()

        for element in self.ep_eff_elements:
            element.calc_eff_stress()
            # element.eff_stress = [-30.e3,-30.e3,0]    #等方圧密
            # element.ep.initial_state_isotropic(element.eff_stress)
            element.ep.initial_state(element.eff_stress)
            element.material.rmu,element.material.rlambda = element.ep.elastic_modulus()
            element.clear_strain()
            element.stress = np.copy(element.eff_stress)

        self._set_ep_initial_state_node_clear()

        for node in self.nodes:
            node.u[:] = np.zeros(self.dof,dtype=np.float64)

        self._set_initial_matrix()

        for element in self.ep_elements:
            element.ep_init_calc_stress_all()
            print("e is",element.ep.e)
        for element in self.ep_eff_elements:
            element.ep_eff_init_calc_stress_all()

    # ---------------------------------------
    def _set_ep_initial_state_node_clear(self):
        for node in self.nodes:
            node.mass = np.zeros(self.dof,dtype=np.float64)
            node.c    = np.zeros(self.dof,dtype=np.float64)
            node.k      = np.zeros(self.dof,dtype=np.float64)
            node.force = np.zeros(self.dof,dtype=np.float64)
            node.static_force = np.zeros(self.dof,dtype=np.float64)
            node.dynamic_force = np.zeros(self.dof,dtype=np.float64)

    # ======================================================================= #
    def update_init(self,dt):
        for node in self.nodes:
            node.inv_mc = 1.0 / (node.mass[:] + 0.5*dt*node.c[:])
            node.mass_inv_mc = node.mass[:]*node.inv_mc[:]
            node.c_inv_mc = node.c[:]*node.inv_mc[:]*0.5*dt
            node.dtdt_inv_mc = dt*dt*node.inv_mc[:]

        self.dt = dt
        self.inv_dt2 = 1./(2.*dt)
        self.inv_dtdt = 1./(dt*dt)

    # ======================================================================= #
    def update_matrix(self):
        for node in self.nodes:
            self._update_matrix_node_init(node)
        for element in self.elements:
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

    def _update_matrix_set_nodes(self,node):
        node.inv_mc = 1.0 / (node.mass[:] + 0.5*self.dt*node.c[:])
        node.mass_inv_mc = node.mass[:]*node.inv_mc[:]
        node.c_inv_mc = node.c[:]*node.inv_mc[:]*0.5*self.dt
        node.dtdt_inv_mc = self.dt*self.dt*node.inv_mc[:]

    # ======================================================================= #
    def update_time(self,acc0,vel0=None,input_wave=False,self_gravity=False,FD=False):
        if FD:
            self.update_matrix()
        else:
            if self_gravity:
                for node in self.nodes:
                    node.dynamic_force = np.copy(node.static_force)
            else:
                for node in self.nodes:
                    node.dynamic_force = np.zeros(self.dof,dtype=np.float64)

        for node in self.nodes:
            self._update_time_node_init(node)

        if input_wave:
            for element in self.input_elements:
                self._update_time_input_wave(element,vel0)
        else:
            for element in self.elements:
                self._update_bodyforce(element,acc0)

        if FD:
            for element in self.e_elements:
                element.mk_B_stress()
                element.mk_cv()
            for element in self.ep_elements:
                element.mk_ep_FD_B_stress()
            for element in self.ep_eff_elements:
                element.mk_ep_FD_eff_B_stress()
        else:
            for element in self.e_elements:
                element.mk_ku_cv()
            for element in self.ep_elements:
                element.mk_ep_B_stress()
            for element in self.ep_eff_elements:
                element.mk_ep_eff_B_stress()

        self._update_time_set_nodes_all()

        for element in self.connected_elements:
            self._update_time_set_connected_elements_(element)
        for element in self.spring_elements:
            self._update_time_set_spring_elements_(element)

        if FD:
            for element in self.output_e_elements:
                element.calc_FD_stress()
        else:
            for element in self.output_e_elements:
                element.calc_stress()

    # ======================================================================= #
    def update_time_input(self,vel0):
        for node in self.node_set:
            node.dynamic_force = np.zeros(self.dof,dtype=np.float64)
            node.force = np.zeros(self.dof,dtype=np.float64)

        for element in self.input_element_set:
            self._update_time_input_wave(element,vel0)

        self._update_mk_internal_force()
        self._update_time_set_nodes_all()

        for element in self.connected_element_set:
            self._update_time_set_connected_elements_(element)
        for element in self.spring_elements:
            self._update_time_set_spring_elements_(element)

        for element in self.output_e_element_set:
            element.calc_stress()
        for element in self.output_ep_element_set:
            element.calc_ep_stress()

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

    def _update_mk_internal_force(self):
        for element in self.elements:
            element.mk_ku_cv()

    def _update_time_set_nodes_all(self):
        for node in self.free_nodes:
            self._update_time_set_free_nodes(node)
        for node in self.fixed_nodes:
            self._update_time_set_fixed_nodes(node)

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
                node.u[i]  = 0.0
            else:
                node.u[i] = node.mass_inv_mc[i]*(2.*u[i]-node.um[i]) + node.c_inv_mc[i]*node.um[i] - node.dtdt_inv_mc[i]*node.force[i]
        node.v[:] = (node.u - node.um) * self.inv_dt2
        node.a[:] = (node.u - 2.*u + node.um) * self.inv_dtdt
        node.um = u

    def _update_time_set_connected_elements_(self,element):
        u = np.zeros_like(element.nodes[0].u)
        a = np.zeros_like(element.nodes[0].a)
        for node in element.nodes:
            u[:] += node.u[:]
            a[:] += node.a[:]
        for node in element.nodes:
            node.u[:] = u[:]/element.nnode
            node.a[:] = a[:]/element.nnode

    def _update_time_set_spring_elements_(self,element):
        u = element.R @ np.hstack(element.u)
        element.f[0] = element.material.kv*(u[0]-u[2])
        element.f[1] = element.material.kh*(u[1]-u[3])

    # ---------------------------------------
    def _multi_update_mk_internal_force(self):
        def mk_ku_cv(element):
            element.mk_ku_cv()
        with concurrent.futures.ProcessPoolExecutor() as executor:
            executor.map(mk_ku_cv,self.elements)

    def _multi_update_time_set_nodes_all(self):
        with concurrent.futures.ProcessPoolExecutor() as executor:
            executor.map(self._update_time_set_free_nodes,self.free_nodes)
        for node in self.fixed_node_set:
            self._update_time_set_fixed_nodes(node)

    # ======================================================================= #
    def print_all(self):
        for node in self.nodes:
            node.print()
        for element in self.elements:
            element.print()
