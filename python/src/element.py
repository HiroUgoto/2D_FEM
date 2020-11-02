import numpy as np

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

    def set_material(self,material):
        self.material = material
        self.rho = material.rho

    # ---------------------------------------------------------
    def set_pointer_list(self):
        self.u, self.v = (), ()
        for node in self.nodes:
            self.u += (node.u.view(),)
            self.v += (node.v.view(),)

    def set_xn(self):       #9*2ノード座標
        self.xn = np.empty([self.nnode,2],dtype=np.float64)
        for i in range(self.nnode):
            self.xn[i,0] = self.nodes[i].xyz[0] + self.u[i][0] # mesh update
            self.xn[i,1] = self.nodes[i].xyz[1] + self.u[i][1] # mesh update
            # self.xn[i,0] = self.nodes[i].xyz[0]
            # self.xn[i,1] = self.nodes[i].xyz[1]

    # ---------------------------------------------------------
    # ---------------------------------------------------------
    def mk_local_matrix_init(self,dof):
        self.dof = dof

        if self.dim == 2:
            self.Mxz = np.empty([self.ng,self.ng,self.dof*self.nnode,self.dof*self.nnode],dtype=np.float64)
            self.Nxz = np.empty([self.ng,self.ng,self.dof,self.dof*self.nnode],dtype=np.float64)
            self.dnxz = self.estyle.init_dn(self.ng)
            self.wxz = np.empty([self.ng,self.ng],dtype=np.float64)

            V = 0.0
            for i,(xi,wx) in enumerate(zip(self.xi,self.w)):        #index付で(xi,w)を取得
                for j,(zeta,wz) in enumerate(zip(self.xi,self.w)):
                    self.Nxz[i,j,:,:] = mk_n(self.dof,self.estyle,self.nnode,xi,zeta)
                    self.Mxz[i,j,:,:] = mk_m(self.Nxz[i,j,:,:])
                    self.dnxz[i,j,:,:] = self.estyle.shape_function_dn(xi,zeta)
                    self.wxz[i,j] = wx*wz

                    det,_ = mk_jacobi(self.estyle,self.xn,self.dnxz[i,j,:,:])
                    detJ = wx*wz*det
                    V += detJ

            self.mass = self.rho*V

        elif self.dim == 1:
            self.Nxz = np.empty([self.ng,self.dof,self.dof*self.nnode],dtype=np.float64)
            self.dnxz = self.estyle.init_dn(self.ng)
            self.wxz = np.empty(self.ng,dtype=np.float64)

            self.imp = self.material.mk_imp(self.dof)
            for i,(xi,wx) in enumerate(zip(self.xi,self.w)):
                self.dnxz[i,:] = self.estyle.shape_function_dn(xi,0.0)
                self.Nxz[i,:,:] = mk_n(self.dof,self.estyle,self.nnode,xi,0.0)
                self.wxz[i] = wx

    # ---------------------------------------------------------
    def mk_local_matrix(self):
        self.M = np.zeros([self.dof*self.nnode,self.dof*self.nnode],dtype=np.float64)
        self.C = np.zeros([self.dof*self.nnode,self.dof*self.nnode],dtype=np.float64)
        self.K = np.zeros([self.dof*self.nnode,self.dof*self.nnode],dtype=np.float64)

        self.De = self.material.mk_d(self.dof)

        self.M_diag = np.diag(self.M)       #対角成分
        self.C_diag = np.diag(self.C)
        self.C_off_diag = np.zeros([self.dof*self.nnode,self.dof*self.nnode],dtype=np.float64)

        if self.dim == 2:
            for i in range(self.ng):
                for j in range(self.ng):
                    det,_ = mk_jacobi(self.estyle,self.xn,self.dnxz[i,j,:,:])
                    detJ = self.wxz[i,j]*det

                    B = mk_b(self.dof,self.estyle,self.nnode,self.xn,self.dnxz[i,j,:,:])
                    K = mk_k(B,self.De)

                    self.M += self.Mxz[i,j,:,:] * detJ
                    self.K += K*detJ

            tr_M = np.trace(self.M)/self.dof
            self.M_diag = np.diag(self.M) * self.mass/tr_M

        elif self.dim == 1:
            if "input" in self.style:
                for i in range(self.ng):
                    det,q = mk_q(self.dof,self.estyle,self.xn,self.dnxz[i,:])
                    detJ = self.wxz[i]*det

                    NqN = mk_nqn(self.dof,self.Nxz[i,:,:],q,self.imp)
                    self.C += NqN*detJ

                self.C_diag = np.diag(self.C)
                self.C_off_diag = self.C - np.diag(self.C_diag)

    # ---------------------------------------------------------
    def mk_local_vector(self):
        self.force = np.zeros(self.dof*self.nnode,dtype=np.float64)

        if self.dof == 1:
            return
        if self.dim == 2:
            V = 0.0
            for i in range(self.ng):
                for j in range(self.ng):
                    det,_ = mk_jacobi(self.estyle,self.xn,self.dnxz[i,j,:,:])
                    detJ = self.wxz[i,j]*det
                    V += detJ
                    self.force += self.Nxz[i,j,1,:]*detJ * self.gravity

            self.force = self.force * self.mass/V

    # ---------------------------------------------------------
    def mk_ku(self):
        ku = np.dot(self.K,np.hstack(self.u))       #横に結合
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += ku[i0:i0+self.dof]

    def mk_cv(self):
        cv = np.dot(self.C_off_diag,np.hstack(self.v))
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += cv[i0:i0+self.dof]

    def mk_ku_cv(self):     #ku項とcv項まとめる
        f = np.dot(self.K,np.hstack(self.u)) + np.dot(self.C_off_diag,np.hstack(self.v))
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += f[i0:i0+self.dof]

    # ---------------------------------------------------------
    def calc_stress(self):
        dn = self.estyle.shape_function_dn(0.0,0.0)
        B = mk_b(self.dof,self.estyle,self.nnode,self.xn,dn)
        self.strain = np.dot(B,np.hstack(self.u))
        self.stress = np.dot(self.De,self.strain)

# ---------------------------------------------------------
def mk_m(N):
    return np.dot(N.T,N)        #18*18

def mk_n(dof,estyle,nnode,xi,zeta):
    n_shape = estyle.shape_function_n(xi,zeta)
    N = np.zeros([dof,dof*nnode],dtype=np.float64)
    e = np.eye(dof)     #単位行列

    for i in range(nnode):
        i0 = dof*i
        N[:,i0:i0+dof] = e*n_shape[i]

    return N

# ---------------------------------------------------------
def mk_nqn(dof,n,q,imp):        #側面境界条件がエネルギー減衰
    qiq = np.dot(np.dot(q.T,imp),q)
    nqn = np.dot(np.dot(n.T,qiq),n)
    return nqn

def mk_q(dof,estyle,xn,dn):
    t = np.dot(xn.T,dn)
    n = np.cross(t,np.array([0.0,0.0,1.0]))
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
    return np.dot(np.dot(B.T,D),B)

def mk_b(dof,estyle,nnode,xn,dn):
    dnj = mk_dnj(estyle,xn,dn)
    if dof == 1:
        B = np.empty([2,nnode],dtype=np.float64)
        for i in range(nnode):
            B[0,i] = dnj[i,0]
            B[1,i] = dnj[i,1]

    elif dof == 2:
        B = np.empty([3,2*nnode],dtype=np.float64)
        for i in range(nnode):
            i0,i1 = 2*i,2*i+1
            B[0,i0],B[0,i1] = dnj[i,0],   0.0
            B[1,i0],B[1,i1] =     0.0 ,dnj[i,1]
            B[2,i0],B[2,i1] = dnj[i,1],dnj[i,0]

    elif dof == 3:
        B = np.zeros([5,3*nnode],dtype=np.float64)
        for i in range(nnode):
            i0,i1,i2 = 3*i,3*i+1,3*i+2
            B[0,i0],B[0,i1] = dnj[i,0],   0.0
            B[1,i0],B[1,i1] =     0.0 ,dnj[i,1]
            B[2,i0],B[2,i1] = dnj[i,1],dnj[i,0]

            B[3,i2] = dnj[i,0]
            B[4,i2] = dnj[i,1]

    return B

# ---------------------------------------------------------
def mk_dnj(estyle,xn,dn):
    _,jacobi_inv = mk_inv_jacobi(estyle,xn,dn)
    return np.dot(dn,jacobi_inv)

def mk_inv_jacobi(estyle,xn,dn):
    det,jacobi = mk_jacobi(estyle,xn,dn)
    jacobi_inv = np.array([ [ jacobi[1,1],-jacobi[0,1]],
                            [-jacobi[1,0], jacobi[0,0]] ]) / det
    return 1.0/det, jacobi_inv

def mk_jacobi(estyle,xn,dn):
    jacobi = np.dot(xn.T,dn)
    det = jacobi[0,0]*jacobi[1,1] - jacobi[0,1]*jacobi[1,0]
    return det,jacobi

# ---------------------------------------------------------
def mk_tau(D,dof,estyle,xn,dn):     #キルヒホッフ応力tau,Kマト算定
    FF,_ = mk_FF(dof,estyle,xn,dn)
    hogehoge = np.empty(dof,dof,dtype=np.float64)

    for i in range(dof):
        for j in range(dof):
            hogehoge[i,j] = 0.5*np.log(FF[i,j])     #Euler対数ひずみ,仮置き
    tau =
    K = np.dot(FF.T,tau)
    return K

def mk_FF(dof,estyle,xn,dn):    #Euler対数ひずみでのBマト
    dnu = mk_dnu(dof,estyle,xn,dn)
    F = np.linalg.inv(np.eye(dof) - dnu)

    FF = np.zeros([2,2],dtype=np.float64)      #FikFjk
    for i in range(dof):
        for j in range(dof):
            FF[i,j] = np.sum(F[i,:]*F[j,:])
    J = np.linalg.det(F)
    return FF,J

def mk_dnu(dof,estyle,dn,u):
    dnu = np.empty(dof,dof,dtype=np.float64)
    dnj = mk_dnj(estyle,xn,dn)
    for i in range(dof):
        for j in range(dof):
            dnu[i,j] = np.sum(u[:][i]*dnj[:,j])
    return dnu
