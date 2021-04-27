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
            M = np.zeros([self.ndof,self.ndof], dtype=np.float64)

            self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)
            self.K = np.zeros([self.ndof,self.ndof],dtype=np.float64)

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
            if "input" in self.style:
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
            if "input" in self.style:
                self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)

                for gp in self.gauss_points:
                    det,q = mk_q(self.dof,self.xnT,gp.dn)
                    detJ = gp.w*det

                    NqN = mk_nqn(self.dof,gp.N,q,self.imp)
                    self.C += NqN*detJ

                self.C_diag = np.diag(self.C)
                self.C_off_diag = self.C - np.diag(self.C_diag)


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
