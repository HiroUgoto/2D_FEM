import numpy as np
import shape_functions

class Element:
    def __init__(self,id,style,param,inode):
        self.id = id
        self.style = style
        self.inode = inode
        self.ng = 5

        self.param = param
        self.set_param(param)

    def print(self):
        print(self.id,":",self.inode,",",self.param)

    def set_nodes(self,nodes):
        self.nodes = nodes

    def set_param(self,param):
        vs,vp,rho = param
        nu = 0.495
        self.rmu = rho/2*vp**2*(1-2*nu)/(1-nu)
        self.rlambda = rho*nu*vp**2/(1-nu)
        self.rho = rho

    # ---------------------------------------------------------
    def mk_local_matrix(self,dof):
        self.dof = dof

        xn = np.empty([self.style,2],dtype=np.float64)
        for i in range(self.style):
            xn[i,0] = self.nodes[i].xyz[0]
            xn[i,1] = self.nodes[i].xyz[1]

        xi_list,w_list = np.polynomial.legendre.leggauss(self.ng)
        self.M = np.zeros([self.dof*self.style,self.dof*self.style],dtype=np.float64)
        self.K = np.zeros([self.dof*self.style,self.dof*self.style],dtype=np.float64)
        self.De = mk_d(self.dof,self.rmu,self.rlambda)

        V = 0.0
        for i in range(self.ng):
            xi,wx = xi_list[i],w_list[i]
            for j in range(self.ng):
                zeta,wz = xi_list[j],w_list[j]

                det,_ = mk_jacobi(self.style,xn,xi,zeta)
                detJ = wx*wz*det

                N = mk_n(self.dof,self.style,xi,zeta)
                M = mk_m(N)
                V += detJ
                self.M += M*detJ

                B = mk_b(self.dof,self.style,xn,xi,zeta)
                K = mk_k(B,self.De)

                self.K += K*detJ

        mass = self.rho*V
        tr_M = np.trace(self.M)/self.dof
        self.M_diag = np.diag(self.M) * mass/tr_M

    # ---------------------------------------------------------
    def mk_ku(self,dof):
        u = np.zeros([dof*self.style])
        for i in range(self.style):
            i0 = dof*i
            u[i0:i0+dof] = self.nodes[i].u[:]

        ku = np.dot(self.K,u)
        for i in range(self.style):
            i0 = dof*i
            self.nodes[i].force[:] += ku[i0:i0+dof]


# ---------------------------------------------------------
def mk_m(N):
    return np.dot(N.T,N)

def mk_n(dof,style,xi,zeta):
    n_shape = shape_functions.shape_function_n(style,xi,zeta)

    N = np.zeros([dof,dof*style],dtype=np.float64)
    e = np.eye(dof)

    for i in range(style):
        i0 = dof*i
        N[:,i0:i0+dof] = e*n_shape[i]

    return N

# ---------------------------------------------------------
def mk_k(B,D):
    return np.dot(np.dot(B.T,D),B)

def mk_d(dof,rmu,rlambda):
    if dof == 1:
        D = np.zeros([2,2],dtype=np.float64)
        D[0,0] = rmu
        D[1,1] = rmu

    elif dof == 2:
        D = np.zeros([3,3],dtype=np.float64)
        D[0,0],D[0,1] = rlambda + 2*rmu, rlambda
        D[1,0],D[1,1] = rlambda, rlambda + 2*rmu
        D[2,2] = rmu

    elif dof == 3:
        D = np.zeros([5,5],dtype=np.float64)
        D[0,0],D[0,1] = rlambda + 2*rmu, rlambda
        D[1,0],D[1,1] = rlambda, rlambda + 2*rmu
        D[2,2] = rmu
        D[3,3] = rmu
        D[4,4] = rmu

    return D

def mk_b(dof,style,xn,xi,zeta):
    dn = mk_dn(style,xn,xi,zeta)

    if dof == 1:
        B = np.zeros([2,style],dtype=np.float64)
        for i in range(style):
            B[0,i] = dn[i,0]
            B[1,i] = dn[i,1]

    elif dof == 2:
        B = np.zeros([3,2*style],dtype=np.float64)
        for i in range(style):
            i0,i1 = 2*i,2*i+1
            B[0,i0],B[0,i1] = dn[i,0],   0.0
            B[1,i0],B[1,i1] =    0.0 ,dn[i,1]
            B[2,i0],B[2,i1] = dn[i,1],dn[i,0]

    elif dof == 3:
        B = np.zeros([5,3*style],dtype=np.float64)
        for i in range(style):
            i0,i1,i2 = 3*i,3*i+1,3*i+2
            B[0,i0],B[0,i1] = dn[i,0],   0.0
            B[1,i0],B[1,i1] =    0.0 ,dn[i,1]
            B[2,i0],B[2,i1] = dn[i,1],dn[i,0]

            B[3,i2] = dn[i,0]
            B[4,i2] = dn[i,1]

    return B

# ---------------------------------------------------------
def mk_dn(style,xn,xi,zeta):
    _,jacobi_inv = mk_inv_jacobi(style,xn,xi,zeta)
    dn = shape_functions.shape_function_dn(style,xi,zeta)

    return np.dot(dn,jacobi_inv)

def mk_inv_jacobi(style,xn,xi,zeta):
    det,jacobi = mk_jacobi(style,xn,xi,zeta)
    jacobi_inv = np.linalg.inv(jacobi)

    return 1.0/det, jacobi_inv

def mk_jacobi(style,xn,xi,zeta):
    dn = shape_functions.shape_function_dn(style,xi,zeta)

    jacobi = np.dot(xn.T,dn)
    det = jacobi[0,0]*jacobi[1,1] - jacobi[0,1]*jacobi[1,0]

    return det,jacobi
    
