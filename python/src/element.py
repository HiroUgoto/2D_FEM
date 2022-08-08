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
            self.ndof += dof
            # du[0]: dun +tension, -compress
            # du[1]: dus +down stream, -reverse
            # theta: compatible with coulomb
            self.du  = np.zeros(dof, dtype=np.float64)
            self.dum = np.zeros(dof, dtype=np.float64)
            self.dv  = np.zeros(dof, dtype=np.float64)

        self.M_diag = np.zeros(self.ndof, dtype=np.float64)

        self.K = np.zeros([self.ndof,self.ndof],dtype=np.float64)
        self.K_diag = np.zeros(self.ndof, dtype=np.float64)
        self.K_off_diag = np.zeros([self.ndof,self.ndof],dtype=np.float64)

        self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)
        self.C_diag = np.zeros(self.ndof, dtype=np.float64)
        self.C_off_diag = np.zeros([self.ndof,self.ndof], dtype=np.float64)

        self.force = np.zeros(self.ndof,dtype=np.float64)

        if self.dim == 2:
            if "X" not in self.style:
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

            # enrich element
            else:
                self.gauss_points = set()
                V = 0.0
                for xi,wx in zip(self.xi,self.w):
                    for zeta,wz in zip(self.xi,self.w):
                        dn = self.estyle.shape_function_dn(xi,zeta)

                        gp = element_style.Gauss_Points(dn,wx*wz,N=None,M=None)
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
            M = np.zeros([self.ndof,self.ndof], dtype=np.float64)

            self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)
            self.K = np.zeros([self.ndof,self.ndof],dtype=np.float64)

            self.De = self.material.mk_d(self.dof)
            self.Dv = self.material.mk_visco(self.dof)

            if "X" not in self.style:
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

            # enrich element
            else:
                self.crack_xp = np.array([3.0,1.0])
                self.crack_theta = np.pi/3.0

                self.crack_edges,self.crack_edges_xi = set_crack_edge(self.xnT,self.crack_xp,self.crack_theta)


                level_set = set_levelset_function(self.xnT,self.crack_xp,self.crack_theta)
                self.Jp,self.Jm = ck_enrich_function_domain_init(level_set)
                self.R = np.array([[np.sin(self.crack_theta), np.cos(self.crack_theta)],
                                   [np.cos(self.crack_theta),-np.sin(self.crack_theta)]])

                ml = 0.0
                for xi,wx in zip(self.xi,self.w):
                    for zeta,wz in zip(self.xi,self.w):
                        dn = self.estyle.shape_function_dn(xi,zeta)
                        n = self.estyle.shape_function_n(xi,zeta)
                        det,dnj = mk_dnj(self.xnT,dn)

                        g,dgj = mk_enrich_function(self.estyle,self.xnT,dn,level_set,self.Jp,self.Jm,xi,zeta)

                        B = mk_b_enrich(self.dof,self.nnode,dnj,dgj,self.R)
                        N = mk_n_enrich(self.dof,self.nnode,n,g,self.R)

                        Me = mk_m(N)
                        K = mk_k(B,self.De)
                        C = mk_k(B,self.Dv)

                        detJ = wx*wz*det
                        M += Me*detJ
                        self.K += K*detJ
                        self.C += C*detJ

                        ml += self.rho*g*g * detJ

                tr_M = np.trace(M)/self.dof
                self.M_diag = np.diag(M) * self.mass/tr_M
                self.M_diag[-self.dof:] = ml

                self.K_diag = np.diag(self.K)
                self.K_off_diag = self.K - np.diag(self.K_diag)

                self.C_diag = np.diag(self.C)
                self.C_off_diag = self.C - np.diag(self.C_diag)


        elif self.dim == 1:
            if "input" or "visco" in self.style:
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
            if "X" not in self.style:
                self.force = np.zeros(self.ndof,dtype=np.float64)
                V = 0.0
                for gp in self.gauss_points:
                    det,_ = mk_jacobi(self.xnT,gp.dn)
                    detJ = gp.w*det
                    V += detJ
                    self.force += gp.N[1,:]*detJ * self.gravity

                self.force = self.force * self.mass/V

            # enrich element
            else:
                self.force = np.zeros(self.ndof,dtype=np.float64)

                level_set = set_levelset_function(self.xnT,self.crack_xp,self.crack_theta)
                Jp,Jm = ck_enrich_function_domain_init(level_set)
                R = np.array([[np.sin(self.crack_theta), np.cos(self.crack_theta)],
                              [np.cos(self.crack_theta),-np.sin(self.crack_theta)]])

                V = 0.0
                for xi,wx in zip(self.xi,self.w):
                    for zeta,wz in zip(self.xi,self.w):
                        dn = self.estyle.shape_function_dn(xi,zeta)
                        n = self.estyle.shape_function_n(xi,zeta)
                        det,dnj = mk_dnj(self.xnT,dn)

                        g,_ = mk_enrich_function(self.estyle,self.xnT,dn,level_set,Jp,Jm,xi,zeta)
                        N = mk_n_enrich(self.dof,self.nnode,n,g,R)

                        detJ = wx*wz*det
                        V += detJ
                        self.force += N[1,:]*detJ * self.gravity

                self.force = self.force * self.mass/V

    # ---------------------------------------------------------
    def mk_local_update(self):
        if self.dim == 2:
            M = np.zeros([self.ndof,self.ndof], dtype=np.float64)
            self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)
            self.Dv = self.material.mk_visco(self.dof)

            gravity_force = np.array([0.0,self.gravity])

            self.force = np.zeros(self.ndof,dtype=np.float64)
            V = 0.0

            for gp in self.gauss_points:
                det,dnj = mk_dnj(self.xnT,gp.dn)
                B = mk_b(self.dof,self.nnode,dnj)
                C = mk_k(B,self.Dv)

                J,F_inv = mk_F_inv(self.nnode,dnj,self.u)
                IFT = np.eye(2) - F_inv.T / J
                sf = IFT @ gravity_force

                detJ = gp.w*det
                M += gp.M*detJ
                self.C += C*detJ

                V += detJ
                self.force += (gp.N.T @ sf) *detJ

            tr_M = np.trace(M)/self.dof
            self.M_diag = np.diag(M) * self.mass/tr_M

            self.C_diag = np.diag(self.C)
            self.C_off_diag = self.C - np.diag(self.C_diag)

            self.force *= self.mass/V

        elif self.dim == 1:
            if "input" or "visco" in self.style:
                self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)

                for gp in self.gauss_points:
                    det,q = mk_q(self.dof,self.xnT,gp.dn)
                    detJ = gp.w*det

                    NqN = mk_nqn(self.dof,gp.N,q,self.imp)
                    self.C += NqN*detJ

                self.C_diag = np.diag(self.C)
                self.C_off_diag = self.C - np.diag(self.C_diag)


    # ---------------------------------------------------------
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
        u = np.append(np.hstack(self.u),self.du)
        v = np.append(np.hstack(self.v),self.dv)

        f = self.K @ u + self.C_off_diag @ v
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += f[i0:i0+self.dof]

        self.enrich_force = f[-self.dof:]

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
        up,um = [],[]
        for xi in self.crack_edges_xi:
            n = self.estyle.shape_function_n(xi[0],xi[1])
            gp = self.estyle.enrich_function_n( 1.0,self.Jp,self.Jm,xi[0],xi[1])
            gm = self.estyle.enrich_function_n(-1.0,self.Jp,self.Jm,xi[0],xi[1])

            Np = mk_n_enrich(self.dof,self.nnode,n,gp,self.R)
            Nm = mk_n_enrich(self.dof,self.nnode,n,gm,self.R)

            ua = np.append(np.hstack(self.u),self.du)
            up += [Np @ np.hstack(ua)]
            um += [Nm @ np.hstack(ua)]

            # print("up:",up)
            # print("um:",um)

        return up,um


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
def mk_nqn(dof,n,q,imp):        #側面境界条件がエネルギー減衰
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

def mk_enrich_function(estyle,xnT,dn,level_set,Jp,Jm,xi,zeta):
    sign = ck_enrich_function_domain(estyle,level_set,xi,zeta)
    g  = estyle.enrich_function_n(sign,Jp,Jm,xi,zeta)
    dg = estyle.enrich_function_dn(sign,Jp,Jm,xi,zeta)

    _,jacobi_inv = mk_inv_jacobi(xnT,dn)
    dgj = dg @ jacobi_inv

    return g, dgj

def mk_b_enrich(dof,nnode,dnj,dgj,R):
    B = mk_b(dof,nnode,dnj)

    dG = np.zeros([3,2],dtype=np.float64)
    dG[0,0] = dgj[0]
    dG[1,1] = dgj[1]
    dG[2,0],dG[2,1] = dgj[1],dgj[0]

    dGR = dG @ R

    B_enrich = np.zeros([3,2*nnode+2],dtype=np.float64)
    B_enrich[:,0:2*nnode] = B[:,:]
    B_enrich[:,2*nnode:] = dGR[:,:]

    return B_enrich

def mk_n_enrich(dof,nnode,n,g,R):
    N_enrich = np.zeros([dof,dof*nnode+dof],dtype=np.float64)

    for i in range(dof):
        N_enrich[i,i:dof*nnode:dof] = n[:]

    N_enrich[:,dof*nnode:] = g * R

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
