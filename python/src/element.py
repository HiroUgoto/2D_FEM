import numpy as np
import sys

import element_style

class Element:
    def __init__(self,id,style,material_id,inode):
        self.id = id
        self.inode = inode
        self.material_id = material_id
        self.gravity = 9.8

        self.set_style(style)

    def print(self):
        print(self.id,":",self.style,",",self.material_id,",",self.inode)

    def set_style(self,style):
        self.style = style
        self.nnode = len(self.inode)

        self.estyle = element_style.set_style(style)
        self.dim = self.estyle.dim

        self.xi,self.w = self.estyle.gauss      #gauss積分点の座標,重み
        self.ng = len(self.xi)      #積分点数

    def set_nodes(self,nodes):
        self.nodes = nodes
        self.node_set = set(nodes)

    def set_material(self,material):
        if material is None:
            self.material = None
            self.rho = None
        else:
            self.material = material
            self.rho = material.rho

    # ---------------------------------------------------------
    def set_pointer_list(self):
        self.u, self.v = (), ()
        for node in self.nodes:
            self.u += (node.u.view(),)
            self.v += (node.v.view(),)

    def set_xn(self):
        self.xnT  = np.empty([2,self.nnode],dtype=np.float64)
        for i in range(self.nnode):
            self.xnT[0,i] = self.nodes[i].xyz[0] + self.u[i][0] # mesh update
            self.xnT[1,i] = self.nodes[i].xyz[1] + self.u[i][1] # mesh update

    # ---------------------------------------------------------
    # ---------------------------------------------------------
    def mk_local_matrix_init(self,dof):
        self.dof = dof
        self.ndof = dof*self.nnode

        if "X" in self.style:
            self.rupture = False
            self.ncrack = 0
            self.crack_edges = []
            self.crack_edges_dummy = []

        self.M_diag = np.zeros(self.ndof, dtype=np.float64)

        self.K = np.zeros([self.ndof,self.ndof],dtype=np.float64)
        self.K_diag = np.zeros(self.ndof, dtype=np.float64)
        self.K_off_diag = np.zeros([self.ndof,self.ndof],dtype=np.float64)

        self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)
        self.C_diag = np.zeros(self.ndof, dtype=np.float64)
        self.C_off_diag = np.zeros([self.ndof,self.ndof], dtype=np.float64)

        self.force = np.zeros(self.ndof,dtype=np.float64)

        if self.dim == 2:
            self.gauss_points = set()
            V = 0.0
            for xi,wx in zip(self.xi,self.w):
                for zeta,wz in zip(self.xi,self.w):
                    N = mk_n(self.dof,self.estyle,self.nnode,xi,zeta)
                    M = mk_m(N)
                    dn = self.estyle.shape_function_dn(xi,zeta)

                    gp = element_style.Gauss_Points(dn,wx*wz,N,M)
                    self.gauss_points.add(gp)

                    det,_ = mk_jacobi(self.xnT,dn)
                    detJ = wx*wz*det
                    V += detJ

            self.mass = self.rho*V

        elif self.dim == 1:
            self.gauss_points = set()
            self.imp = self.material.mk_imp(self.dof)

            for xi,wx in zip(self.xi,self.w):
                dn = self.estyle.shape_function_dn(xi,0.0)
                N = mk_n(self.dof,self.estyle,self.nnode,xi,0.0)
                w = wx

                gp = element_style.Gauss_Points(dn,wx,N)
                self.gauss_points.add(gp)

    # ---------------------------------------------------------
    def mk_local_matrix(self):
        if self.dim == 2:
            if "X" in self.style and self.rupture:
                self.mk_local_matrix_enrich()
            else:
                M = np.zeros([self.ndof,self.ndof], dtype=np.float64)

                self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)
                self.K = np.zeros([self.ndof,self.ndof], dtype=np.float64)

                self.De = self.material.mk_d(self.dof)
                self.Dv = self.material.mk_visco(self.dof)

                for gp in self.gauss_points:
                    det,dnj = mk_dnj(self.xnT,gp.dn)
                    B = mk_b(self.dof,self.nnode,dnj)
                    K = mk_k(B,self.De)
                    C = mk_k(B,self.Dv)

                    detJ = gp.w*det
                    M += gp.M*detJ
                    self.K += K*detJ
                    self.C += C*detJ

                tr_M = np.trace(M)/self.dof
                self.M_diag = np.diag(M) * self.mass/tr_M

                self.K_diag = np.diag(self.K)
                self.K_off_diag = self.K - np.diag(self.K_diag)

                self.C_diag = np.diag(self.C)
                self.C_off_diag = self.C - np.diag(self.C_diag)

        elif self.dim == 1:
            if ("input" in self.style) or ("visco" in self.style):
                self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)

                for gp in self.gauss_points:
                    det,q = mk_q(self.dof,self.xnT,gp.dn)
                    detJ = gp.w*det

                    NqN = mk_nqn(self.dof,gp.N,q,self.imp)
                    self.C += NqN*detJ

                self.C_diag = np.diag(self.C)
                self.C_off_diag = self.C - np.diag(self.C_diag)

    # ---------------------------------------------------------
    def mk_local_vector(self):
        # if self.dof == 1:
        #     return
        if self.dim == 2:
            if ("X" in self.style) and self.rupture:
                self.mk_local_vector_enrich()
            else:
                self.force = np.zeros(self.ndof,dtype=np.float64)
                V = 0.0
                for gp in self.gauss_points:
                    det,_ = mk_jacobi(self.xnT,gp.dn)
                    detJ = gp.w*det
                    V += detJ
                    self.force += gp.N[1,:]*detJ * self.gravity

                self.force = self.force * self.mass/V

    # ---------------------------------------------------------
    def mk_local_update(self):
        if self.dim == 2:
            if "X" in self.style and self.rupture:
                self.mk_local_matrix_enrich()
                self.mk_local_vector_enrich()
            else:
                M = np.zeros([self.ndof,self.ndof], dtype=np.float64)

                self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)
                self.K = np.zeros([self.ndof,self.ndof],dtype=np.float64)

                self.De = self.material.mk_d(self.dof)
                self.Dv = self.material.mk_visco(self.dof)

                gravity_force = np.array([0.0,self.gravity])
                self.force = np.zeros(self.ndof,dtype=np.float64)
                V = 0.0
                for gp in self.gauss_points:
                    det,dnj = mk_dnj(self.xnT,gp.dn)
                    B = mk_b(self.dof,self.nnode,dnj)
                    K = mk_k(B,self.De)
                    C = mk_k(B,self.Dv)

                    J,F_inv = mk_F_inv(self.nnode,dnj,self.u)
                    IFT = np.eye(2) - F_inv.T / J
                    sf = IFT @ gravity_force

                    detJ = gp.w*det
                    M += gp.M*detJ
                    self.K += K*detJ
                    self.C += C*detJ

                    V += detJ
                    self.force += (gp.N.T @ sf) *detJ

                tr_M = np.trace(M)/self.dof
                self.M_diag = np.diag(M) * self.mass/tr_M

                self.K_diag = np.diag(self.K)
                self.K_off_diag = self.K - np.diag(self.K_diag)

                self.C_diag = np.diag(self.C)
                self.C_off_diag = self.C - np.diag(self.C_diag)

                self.force *= self.mass/V


        elif self.dim == 1:
            if ("input" in self.style) or ("visco" in self.style):
                self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)

                for gp in self.gauss_points:
                    det,q = mk_q(self.dof,self.xnT,gp.dn)
                    detJ = gp.w*det

                    NqN = mk_nqn(self.dof,gp.N,q,self.imp)
                    self.C += NqN*detJ

                self.C_diag = np.diag(self.C)
                self.C_off_diag = self.C - np.diag(self.C_diag)


    # ---------------------------------------------------------
    def check_enrich_rupture(self,all_crack_edges):
        self.ft = 30000.0

        if not self.rupture:
            self.calc_stress()

            # check tension crack
            stress = np.array([ [self.stress[0], self.stress[2]],
                                [self.stress[2], self.stress[1]]])
            s,v = np.linalg.eig(stress)
            s_max = s[0]
            n_max = v[:,0]

            if self.ft < s_max:
                self.ncrack = 2
                self.ndof = self.nnode*self.dof + self.ncrack*self.dof

                # du[0]: dun +tension, -compress
                # du[1]: dus +down stream, -reverse
                # theta: compatible with coulomb
                self.du_list = []
                self.dum_list = []
                self.dv_list = []
                for ic in range(self.ncrack):
                    self.du_list += [np.zeros(self.dof, dtype=np.float64)]
                    self.dum_list += [np.zeros(self.dof, dtype=np.float64)]
                    self.dv_list += [np.zeros(self.dof, dtype=np.float64)]

                theta0 = np.arctan2( n_max[0], n_max[1])
                theta1 = np.arctan2( n_max[0], n_max[1])
                self.crack_theta_list = [theta0,theta1]
                self.crack_open_side = [0.0, np.pi]

                # initial crack
                n = self.estyle.shape_function_n(0.0,0.0)
                xp = self.xnT @ n
                self.crack_xp_list = [xp,xp]
                self.rupture = True
                self.crack_edges,_ = set_crack_edge(self.xnT,xp,theta0)

                # continuous crack
                if all_crack_edges:
                    for c in all_crack_edges:
                        is_inside,xi = self.check_inside(c,margin=0.01)
                        if is_inside:
                            self.ncrack = 0
                            self.ndof = self.nnode*self.dof + self.ncrack*self.dof
                            theta = np.arctan2( n_max[0], n_max[1])

                            self.crack_theta_list = []
                            self.crack_edges = []

                            self.crack_theta_list_dummy = [theta]
                            self.crack_edges_dummy = [c]

                            self.rupture = True
                            return True
                return True

        return False

    # ---------------------------------------------------------
    def make_half_crack(self,all_crack_edges,all_crack_theta):
        if self.ncrack == 0:
            if all_crack_edges:
                for ic in range(len(all_crack_edges)):
                    c = all_crack_edges[ic]
                    theta = all_crack_theta[ic]
                    is_inside,xi = self.check_inside(c,margin=0.01)
                    if is_inside:
                        self.ncrack = 1
                        self.ndof = self.nnode*self.dof + self.ncrack*self.dof

                        self.du_list = [np.zeros(self.dof, dtype=np.float64)]
                        self.dum_list = [np.zeros(self.dof, dtype=np.float64)]
                        self.dv_list = [np.zeros(self.dof, dtype=np.float64)]

                        d,cxi = set_crack_edge(self.xnT,c,theta)
                        xp = 0.5*(d[0]+d[1])
                        self.crack_xp_list = [xp]
                        self.crack_theta_list = [theta]

                        u = np.array([np.cos(theta),-np.sin(theta)])
                        if np.dot(u,(c-xp)) > 0.0:
                            self.crack_open_side = [0.0]
                        else:
                            self.crack_open_side = [np.pi]

                        self.rupture = True
                        self.crack_edges = [ 2*xp-c ]

                        self.crack_theta_list_dummy = []
                        self.crack_edges_dummy = []

    # ---------------------------------------------------------
    def make_connect_crack(self,all_crack_edges_dummy,all_crack_theta_dummy):
        if self.ncrack == 1:
            if all_crack_edges_dummy:
                for ic in range(len(all_crack_edges_dummy)):
                    c = all_crack_edges_dummy[ic]
                    theta = all_crack_theta_dummy[ic]
                    is_inside,xi = self.check_inside(c,margin=0.01)
                    if is_inside:
                        self.ncrack = 2
                        self.ndof = self.nnode*self.dof + self.ncrack*self.dof

                        self.du_list += [np.zeros(self.dof, dtype=np.float64)]
                        self.dum_list += [np.zeros(self.dof, dtype=np.float64)]
                        self.dv_list += [np.zeros(self.dof, dtype=np.float64)]

                        xp = self.crack_xp_list[0]
                        self.crack_xp_list += [xp]
                        self.crack_theta_list += [theta]

                        u = np.array([np.cos(theta),-np.sin(theta)])
                        if np.dot(u,(c-xp)) > 0.0:
                            self.crack_open_side += [0.0]
                        else:
                            self.crack_open_side += [np.pi]

                        d,cxi = set_crack_edge(self.xnT,xp,theta)
                        if np.dot((c-xp),(d[0]-xp)) > 0.0:
                            self.crack_edges = [d[0]]
                        else:
                            self.crack_edges = [d[1]]


    # ---------------------------------------------------------
    def mk_local_matrix_enrich(self):
        M = np.zeros([self.ndof,self.ndof], dtype=np.float64)

        self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)
        self.K = np.zeros([self.ndof,self.ndof], dtype=np.float64)

        self.De = self.material.mk_d(self.dof)
        self.Dv = self.material.mk_visco(self.dof)

        self.level_set_list = []
        self.Jp_list,self.Jm_list = [],[]
        self.R_list = []
        for ic in range(self.ncrack):
            crack_xp = self.crack_xp_list[ic]
            crack_theta = self.crack_theta_list[ic]

            level_set = set_levelset_function(self.xnT,crack_xp,crack_theta)

            Jp,Jm = ck_enrich_function_domain_init(level_set)
            R = np.array([[np.sin(crack_theta), np.cos(crack_theta)],
                          [np.cos(crack_theta),-np.sin(crack_theta)]])

            self.level_set_list += [level_set]
            self.Jp_list += [Jp]
            self.Jm_list += [Jm]
            self.R_list += [R]

        ml = 0.0
        for xi,wx in zip(self.xi,self.w):
            for zeta,wz in zip(self.xi,self.w):
                dn = self.estyle.shape_function_dn(xi,zeta)
                n = self.estyle.shape_function_n(xi,zeta)
                det,dnj = mk_dnj(self.xnT,dn)

                gh_list,dgh_list = [],[]
                sum_gh = 0.0
                for ic in range(self.ncrack):
                    g,dgj = mk_enrich_function(self.estyle,self.xnT,dn,self.level_set_list[ic],self.Jp_list[ic],self.Jm_list[ic],xi,zeta)

                    x = self.xnT @ n
                    h,dhj = mk_tip_enrich_function(x,self.crack_xp_list[ic],self.crack_theta_list[ic],self.crack_open_side[ic])

                    gh = g*h
                    dgh = dgj*h + dhj*g

                    gh_list += [gh]
                    dgh_list += [dgh]
                    sum_gh += gh


                B = mk_b_enrich(self.dof,self.nnode,dnj,self.ncrack,dgh_list,self.R_list)
                N = mk_n_enrich(self.dof,self.nnode,n,self.ncrack,gh_list,self.R_list)

                Me = mk_m(N)
                K = mk_k(B,self.De)
                C = mk_k(B,self.Dv)

                detJ = wx*wz*det
                M += Me*detJ
                self.K += K*detJ
                self.C += C*detJ

                ml += self.rho*sum_gh**2 * detJ

        tr_M = np.trace(M[:self.dof*self.nnode,:self.dof*self.nnode])/self.dof
        self.M_diag = np.diag(M) * self.mass/tr_M

        self.M_diag[self.dof*self.nnode:self.dof*self.nnode+self.dof*self.ncrack] = ml/self.ncrack

        self.K_diag = np.diag(self.K)
        self.K_off_diag = self.K - np.diag(self.K_diag)

        self.C_diag = np.diag(self.C)
        self.C_off_diag = self.C - np.diag(self.C_diag)

    # ---------------------------------------------------------
    def mk_local_vector_enrich(self):
        self.force = np.zeros(self.ndof,dtype=np.float64)

        V = 0.0
        for xi,wx in zip(self.xi,self.w):
            for zeta,wz in zip(self.xi,self.w):
                dn = self.estyle.shape_function_dn(xi,zeta)
                n = self.estyle.shape_function_n(xi,zeta)
                det,dnj = mk_dnj(self.xnT,dn)

                gh_list = []
                for ic in range(self.ncrack):
                    g,_ = mk_enrich_function(self.estyle,self.xnT,dn,self.level_set_list[ic],self.Jp_list[ic],self.Jm_list[ic],xi,zeta)

                    x = self.xnT @ n
                    h,_ = mk_tip_enrich_function(x,self.crack_xp_list[ic],self.crack_theta_list[ic],self.crack_open_side[ic])

                    gh = g*h
                    gh_list += [gh]

                N = mk_n_enrich(self.dof,self.nnode,n,self.ncrack,gh_list,self.R_list)

                detJ = wx*wz*det
                V += detJ
                self.force += N[1,:]*detJ * self.gravity

        self.force = self.force * self.mass/V

    # ---------------------------------------------------------



    # ---------------------------------------------------------
    def mk_ku(self):
        ku = self.K @ np.hstack(self.u)
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += ku[i0:i0+self.dof]

    def mk_ku_u(self,u):
        ku = self.K @ np.hstack(u)
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += ku[i0:i0+self.dof]

    def mk_cv(self):
        cv = self.C_off_diag @ np.hstack(self.v)
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += cv[i0:i0+self.dof]

    def mk_ku_cv(self):
        f = self.K @ np.hstack(self.u) + self.C_off_diag @ np.hstack(self.v)
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += f[i0:i0+self.dof]

    def mk_bodyforce(self,acc0):
        if self.dim == 2:
            self.force = np.zeros(self.ndof,dtype=np.float64)
            V = 0.0
            for gp in self.gauss_points:
                det,_ = mk_jacobi(self.xnT,gp.dn)
                detJ = gp.w*det
                V += detJ
                self.force += (gp.N[0,:]*acc0[0] + gp.N[1,:]*acc0[1])*detJ

            self.force = self.force * self.mass/V

    # -------------------------------------------------------
    def mk_ku_cv_enrich(self):
        if self.rupture:
            u = np.append(np.hstack(self.u),self.du_list)
            v = np.append(np.hstack(self.v),self.dv_list)

            f = self.K @ u + self.C_off_diag @ v
            for i in range(self.nnode):
                i0 = self.dof*i
                self.nodes[i].force[:] += f[i0:i0+self.dof]

            self.enrich_force = f[-self.dof*self.ncrack:]

        else:
            self.mk_ku_cv()

    # --------------------------------------------------------
    def mk_B_stress(self):
        if self.dim == 1:
            self.mk_ku()

        elif self.dim == 2:
            force = np.zeros(self.ndof,dtype=np.float64)

            for gp in self.gauss_points:
                det,dnj = mk_dnj(self.xnT,gp.dn)
                BT = mk_b_T(self.dof,self.nnode,dnj)
                stress = Hencky_stress(self.dof,self.nnode,self.De,dnj,self.u)

                detJ = gp.w*det
                force += BT @ stress * detJ

            for i in range(self.nnode):
                i0 = self.dof*i
                self.nodes[i].force[:] += force[i0:i0+self.dof]

    def mk_B_stress_u(self,u):
        if self.dim == 1:
            self.mk_ku_u(u)

        elif self.dim == 2:
            force = np.zeros(self.ndof,dtype=np.float64)

            for gp in self.gauss_points:
                det,dnj = mk_dnj(self.xnT,gp.dn)
                BT = mk_b_T(self.dof,self.nnode,dnj)
                stress = Hencky_stress(self.dof,self.nnode,self.De,dnj,u)

                detJ = gp.w*det
                force += BT @ stress * detJ

            for i in range(self.nnode):
                i0 = self.dof*i
                self.nodes[i].force[:] += force[i0:i0+self.dof]


    # ---------------------------------------------------------
    def calc_stress(self):
        dn = self.estyle.shape_function_dn(0.0,0.0)
        _,dnj = mk_dnj(self.xnT,dn)
        B = mk_b(self.dof,self.nnode,dnj)
        self.strain = B @ np.hstack(self.u)
        self.stress = self.De @ self.strain


    # ---------------------------------------------------------
    def calc_crack_edge_disp(self):
        crack_list = []
        for ic in range(self.ncrack):
            crack_dict = {}
            crack_dict["x"] = []
            crack_dict["up"] = []
            crack_dict["um"] = []

            crack_xp = self.crack_xp_list[ic]
            _,crack_xip = self.check_inside(crack_xp)
            crack_theta = self.crack_theta_list[ic]
            crack_edges,crack_edges_xi = set_crack_edge(self.xnT,crack_xp,crack_theta)

            n = self.estyle.shape_function_n(crack_xip[0],crack_xip[1])
            x = self.xnT @ n

            ghp_list,ghm_list = [],[]
            for jc in range(self.ncrack):
                gp = self.estyle.enrich_function_n( 1.0,self.Jp_list[jc],self.Jm_list[jc],crack_xip[0],crack_xip[1])
                gm = self.estyle.enrich_function_n(-1.0,self.Jp_list[jc],self.Jm_list[jc],crack_xip[0],crack_xip[1])

                ghp_list += [gp*0.5]
                ghm_list += [gm*0.5]

            Np = mk_n_enrich(self.dof,self.nnode,n,self.ncrack,ghp_list,self.R_list)
            Nm = mk_n_enrich(self.dof,self.nnode,n,self.ncrack,ghm_list,self.R_list)

            ua = np.append(np.hstack(self.u),self.du_list)
            up = Np @ np.hstack(ua)
            um = Nm @ np.hstack(ua)

            crack_dict["x"] += [x]
            crack_dict["up"] += [up]
            crack_dict["um"] += [um]

            for xi in crack_edges_xi:
                n = self.estyle.shape_function_n(xi[0],xi[1])
                x = self.xnT @ n

                h = mk_tip_enrich_function_h(x,crack_xp,crack_theta,self.crack_open_side[ic])

                if h > 0.95:
                    ghp_list,ghm_list = [],[]
                    for jc in range(self.ncrack):
                        gp = self.estyle.enrich_function_n( 1.0,self.Jp_list[jc],self.Jm_list[jc],xi[0],xi[1])
                        gm = self.estyle.enrich_function_n(-1.0,self.Jp_list[jc],self.Jm_list[jc],xi[0],xi[1])
                        hc = mk_tip_enrich_function_h(x,self.crack_xp_list[jc],self.crack_theta_list[jc],self.crack_open_side[jc])

                        ghp_list += [gp*hc]
                        ghm_list += [gm*hc]

                    Np = mk_n_enrich(self.dof,self.nnode,n,self.ncrack,ghp_list,self.R_list)
                    Nm = mk_n_enrich(self.dof,self.nnode,n,self.ncrack,ghm_list,self.R_list)

                    ua = np.append(np.hstack(self.u),self.du_list)
                    up = Np @ np.hstack(ua)
                    um = Nm @ np.hstack(ua)

                    crack_dict["x"] += [x]
                    crack_dict["up"] += [up]
                    crack_dict["um"] += [um]

            crack_list += [crack_dict]

        return crack_list

    # ---------------------------------------------------------
    def check_inside(self,x,margin=0.0):
        def J_func(xi,x,xnT,shape_function_n,shape_function_dn):
            n = shape_function_n(xi[0],xi[1])
            return xnT@n - x

        def dJ_func(xi,x,xnT,shape_function_n,shape_function_dn):
            dn = shape_function_dn(xi[0],xi[1])
            _,dJ = mk_jacobi(xnT,dn)
            return dJ

        xi = np.zeros(2)
        for itr in range(20):
            n = self.estyle.shape_function_n(xi[0],xi[1])
            dn = self.estyle.shape_function_dn(xi[0],xi[1])

            J = self.xnT@n - x
            _,dJ = mk_jacobi(self.xnT,dn)

            r = np.linalg.solve(dJ,J)
            if np.linalg.norm(r) < 1e-8:
                break

            xi -= r

        if (-1.0-margin <= xi[0] <= 1.0+margin) and (-1.0-margin <= xi[1] <= 1.0+margin):
            is_inside = True
        else:
            is_inside = False

        return is_inside,xi

# ---------------------------------------------------------
def mk_m(N):
    return N.T @ N

def mk_n(dof,estyle,nnode,xi,zeta):
    n_shape = estyle.shape_function_n(xi,zeta)
    N = np.zeros([dof,dof*nnode],dtype=np.float64)

    for i in range(dof):
        N[i,i::dof] = n_shape[:]

    return N

# ---------------------------------------------------------
def mk_nqn(dof,n,q,imp):
    nqn = np.linalg.multi_dot([n.T,q.T,imp,q,n])
    return nqn

def mk_q(dof,xnT,dn):
    t = xnT @ dn
    n = np.cross(t,[0.0,0.0,1.0])
    det = np.linalg.norm(t)

    if dof == 1:
        q = np.array([1.0])
    elif dof == 2:
        q = np.array([[n[0],n[1]],
                      [t[0],t[1]]]) / det
    elif dof == 3:
        q = np.array([[n[0],n[1],0.0],
                      [t[0],t[1],0.0],
                      [0.0 ,0.0 ,1.0]]) / det

    return det, q

# ---------------------------------------------------------
def mk_k(B,D):
    return np.linalg.multi_dot([B.T,D,B])

def mk_b(dof,nnode,dnj):
    if dof == 1:
        B = np.zeros([2,nnode],dtype=np.float64)
        B[0,:] = dnj[:,0]
        B[1,:] = dnj[:,1]

    elif dof == 2:
        B = np.zeros([3,2*nnode],dtype=np.float64)
        B[0,0::2] = dnj[:,0]
        B[1,1::2] = dnj[:,1]
        B[2,0::2],B[2,1::2] = dnj[:,1],dnj[:,0]

    elif dof == 3:
        B = np.zeros([5,3*nnode],dtype=np.float64)
        B[0,0::3] = dnj[:,0]
        B[1,1::3] = dnj[:,1]
        B[2,0::3],B[2,1::3] = dnj[:,1],dnj[:,0]

        B[3,2::3] = dnj[:,0]
        B[4,2::3] = dnj[:,1]

    return B

def mk_b_T(dof,nnode,dnj):
    if dof == 1:
        B = np.zeros([nnode,2],dtype=np.float64)
        B[:,0] = dnj[:,0]
        B[:,1] = dnj[:,1]

    elif dof == 2:
        B = np.zeros([2*nnode,3],dtype=np.float64)
        B[0::2,0] = dnj[:,0]
        B[1::2,1] = dnj[:,1]
        B[0::2,2],B[1::2,2] = dnj[:,1],dnj[:,0]

    elif dof == 3:
        B = np.zeros([3*nnode,5],dtype=np.float64)
        B[0::3,0] = dnj[:,0]
        B[1::3,1] = dnj[:,1]
        B[0::3,2],B[1::3,2] = dnj[:,1],dnj[:,0]

        B[2::3,3] = dnj[:,0]
        B[2::3,4] = dnj[:,1]

    return B

# ---------------------------------------------------------
def set_levelset_function(xnT,xp,theta):
    u = np.array([np.cos(theta),-np.sin(theta)])
    # if -np.sin(theta) > 0.0:
    #     u = -u

    level_set = []
    for x in xnT.T:
        px = x - xp
        level_set += [np.cross(u,px)]

    return np.array(level_set)

def set_crack_edge(xnT,xp,theta):
    level_set = set_levelset_function(xnT,xp,theta)

    crack_edges = []
    crack_edges_xi = []
    for i in range(len(level_set)):
        i0 = i
        if i+1 < len(level_set):
            i1 = i+1
        else:
            i1 = 0

        l0 = level_set[i0]
        l1 = level_set[i1]

        if l0*l1 < 0.0:
            t = (l0+l1)/(l0-l1)
            xe = 0.5*(xnT[:,i1]-xnT[:,i0])*t + 0.5*(xnT[:,i1]+xnT[:,i0])
            crack_edges += [xe]

            if i == 0:
                xi,zeta = t,-1.0
            elif i == 1:
                xi,zeta = 1.0,t
            elif i == 2:
                xi,zeta = -t,1.0
            elif i == 3:
                xi,zeta = -1.0,-t

            crack_edges_xi += [np.array([xi,zeta])]

    return crack_edges, crack_edges_xi

def ck_enrich_function_domain_init(level_set):
    Jp,Jm = [],[]

    for id in range(len(level_set)):
        if level_set[id] > 0:
            Jp += [id]
        else:
            Jm += [id]

    return Jp,Jm

def ck_enrich_function_domain(estyle,level_set,xi,zeta):
    n_shape = estyle.shape_function_n(xi,zeta)
    sign = level_set @ n_shape

    return sign


def mk_tip_enrich_function_h(x,xp,theta,open_side):
    x_xp = x-xp
    r = np.linalg.norm(x_xp)
    u = np.array([np.cos(theta+open_side),-np.sin(theta+open_side)])
    p = np.dot(x_xp,u)/r

    if p >= 1.0:
        phi = 0.0
    elif p <= -1.0:
        phi = np.pi
    else:
        phi = np.arccos(p)
    h = 1.0 - phi/np.pi

    return h

def mk_tip_enrich_function(x,xp,theta,open_side):
    x_xp = x-xp
    r = np.linalg.norm(x_xp)
    u = np.array([np.cos(theta+open_side),-np.sin(theta+open_side)])
    p = np.dot(x_xp,u)/r

    p_x = u[0]/r - x_xp[0]*np.dot(x_xp,u)/r**3
    p_z = u[1]/r - x_xp[1]*np.dot(x_xp,u)/r**3

    phi = np.arccos(p)
    h = 1.0 - phi/np.pi
    h_x = p_x / np.pi / np.sqrt(1-p**2)
    h_z = p_z / np.pi / np.sqrt(1-p**2)

    return h, np.array([h_x,h_z])


def mk_enrich_function(estyle,xnT,dn,level_set,Jp,Jm,xi,zeta):
    sign = ck_enrich_function_domain(estyle,level_set,xi,zeta)
    g  = estyle.enrich_function_n(sign,Jp,Jm,xi,zeta)
    dg = estyle.enrich_function_dn(sign,Jp,Jm,xi,zeta)

    _,jacobi_inv = mk_inv_jacobi(xnT,dn)
    dgj = dg @ jacobi_inv

    return g, dgj

def mk_b_enrich(dof,nnode,dnj,ncrack,dgj,R):
    B = mk_b(dof,nnode,dnj)

    dGR_list = []
    for ic in range(ncrack):
        dG = np.zeros([3,2],dtype=np.float64)
        dG[0,0] = dgj[ic][0]
        dG[1,1] = dgj[ic][1]
        dG[2,0],dG[2,1] = dgj[ic][1],dgj[ic][0]
        dGR_list += [dG @ R[ic]]

    B_enrich = np.zeros([3,2*nnode+2*ncrack],dtype=np.float64)
    B_enrich[:,0:2*nnode] = B[:,:]

    for ic in range(ncrack):
        B_enrich[:,2*nnode+2*ic:2*nnode+2*(ic+1)] = dGR_list[ic][:,:]

    return B_enrich

def mk_n_enrich(dof,nnode,n,ncrack,g,R):
    N_enrich = np.zeros([dof,dof*nnode+dof*ncrack],dtype=np.float64)

    for i in range(dof):
        N_enrich[i,i:dof*nnode:dof] = n[:]

    for ic in range(ncrack):
        N_enrich[:,dof*nnode+dof*ic:dof*nnode+dof*(ic+1)] = g[ic] * R[ic]

    return N_enrich

# ---------------------------------------------------------
def Hencky_stress(dof,nnode,D,dnj,u):
    J,strain = Euler_log_strain(nnode,dnj,u)
    # strain = micro_strain(nnode,dnj,u)
    # J,_ = mk_F(nnode,xnT,dn,u)

    strain_vector = [strain[0,0],strain[1,1],strain[0,1]+strain[1,0]]
    K_stress = np.matmul(D,strain_vector)

    return K_stress/J

def Euler_log_strain(nnode,dnj,u):
    J,FF = mk_FF(nnode,dnj,u)
    L,P = np.linalg.eigh(FF)
    log_L = np.log(L)

    EL_strain = 0.5* P @ np.diag(log_L) @ P.T
    return J,EL_strain

def micro_strain(nnode,dnj,u):
    dnu = mk_dnu(nnode,dnj,u)
    strain = np.array([ [dnu[0,0], 0.5*(dnu[0,1]+dnu[1,0])],
                        [0.5*(dnu[0,1]+dnu[1,0]), dnu[1,1]] ])
    return strain

def mk_FF(nnode,dnj,u):
    J,F = mk_F(nnode,dnj,u)
    FF = np.matmul(F,F.T)
    return J,FF

def mk_F_inv(nnode,dnj,u):
    dnu = mk_dnu(nnode,dnj,u)
    F_inv = np.array([ [1.0-dnu[0,0],     -dnu[0,1]],
                     [     -dnu[1,0],  1.0-dnu[1,1]] ])
    det = (1.0-dnu[0,0])*(1.0-dnu[1,1]) - dnu[0,1]*dnu[1,0]
    return 1./det, F_inv

def mk_F(nnode,dnj,u):
    dnu = mk_dnu(nnode,dnj,u)
    det = (1.0-dnu[0,0])*(1.0-dnu[1,1]) - dnu[0,1]*dnu[1,0]
    F = np.array([ [1.0-dnu[1,1],     dnu[0,1]],
                 [    dnu[1,0], 1.0-dnu[0,0]] ]) / det
    return 1./det, F

def mk_dnu(nnode,dnj,u):
    u_mt = np.array(u)
    return np.matmul(u_mt.T,dnj)

# ---------------------------------------------------------
def mk_dnj(xnT,dn):
    det,jacobi_inv = mk_inv_jacobi(xnT,dn)
    return det, np.matmul(dn,jacobi_inv)

def mk_inv_jacobi(xnT,dn):
    det,jacobi = mk_jacobi(xnT,dn)
    jacobi_inv = np.array([ [ jacobi[1,1],-jacobi[0,1]],
                            [-jacobi[1,0], jacobi[0,0]] ]) / det
    return det, jacobi_inv

def mk_jacobi(xnT,dn):
    jacobi = np.matmul(xnT,dn)
    det = jacobi[0,0]*jacobi[1,1] - jacobi[0,1]*jacobi[1,0]
    return det, jacobi
