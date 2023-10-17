import numpy as np

class Material:
    def __init__(self,id,style,param):
        self.id = id
        self.style = style

        self.set_param(param)

    def print(self):
        print(self.id,":",self.style,self.param)

    def set_param(self,param):
        self.param = param

        if self.style == "vs_vp_rho":
            vs,vp,rho = param
            self.rmu = rho*vs*vs
            self.rlambda = rho*vp*vp - 2.*self.rmu
            self.rho = rho

        elif self.style == "nu_vp_rho":
            nu,vp,rho = param
            self.rmu = rho/2*vp**2*(1-2*nu)/(1-nu)
            self.rlambda = rho*nu*vp**2/(1-nu)
            self.rho = rho

        elif self.style == "nu_vs_rho":
            nu,vs,rho = param
            self.rmu = rho*vs**2
            self.rlambda = 2*nu/(1-2*nu) * self.rmu
            self.rho = rho

        elif self.style == "nu_E_rho":
            nu,E,rho = param
            self.rmu = E/2.0/(1+nu)
            self.rlambda = 2*nu/(1-2*nu) * self.rmu
            self.rho = rho

        elif self.style == "spring_normal":
            self.rho = 0.0
            self.kv,self.kh,n0,n1 = param

            norm = np.sqrt(n0*n0+n1*n1)
            n0,n1 = n0/norm,n1/norm

            self.R = np.zeros([2,2], dtype=np.float64)
            self.R[0,0],self.R[0,1] =  n0, n1
            self.R[1,0],self.R[1,1] = -n1, n0

        elif self.style == "slider_normal":
            self.rho = 0.0
            n0,n1 = param

            norm = np.sqrt(n0*n0+n1*n1)
            n0,n1 = n0/norm,n1/norm

            self.R = np.zeros([2,2], dtype=np.float64)
            self.R[0,0],self.R[0,1] =  n0, n1
            self.R[1,0],self.R[1,1] = -n1, n0

        elif self.style == "ep_Li":
            self.rho = param[0]
            self.param = param[1:]

        elif self.style == "ep_eff_Li":
            self.rho = param[0]
            self.param = param[1:]

            self.rho_w = 1000.0
            self.Kw = 2.25e9

            e0 = param[4]
            n0 = e0/(1 + e0)
            self.rho = (1-n0)*param[0] + n0*self.rho_w
            # print("rho_t:",self.rho)

        elif self.style == "ep_GHES":
            self.rho = param[0]
            self.param = param[1:]

    # ---------------------------------------------------------
    def mk_d(self,dof):
        if dof == 1:
            D = np.zeros([2,2],dtype=np.float64)
            D[0,0] = self.rmu
            D[1,1] = self.rmu

        elif dof == 2:
            D = np.zeros([3,3],dtype=np.float64)
            D[0,0],D[0,1] = self.rlambda + 2*self.rmu, self.rlambda
            D[1,0],D[1,1] = self.rlambda, self.rlambda + 2*self.rmu
            D[2,2] = self.rmu

        elif dof == 3:
            D = np.zeros([5,5],dtype=np.float64)
            D[0,0],D[0,1] = self.rlambda + 2*self.rmu, self.rlambda
            D[1,0],D[1,1] = self.rlambda, self.rlambda + 2*self.rmu
            D[2,2] = self.rmu
            D[3,3] = self.rmu
            D[4,4] = self.rmu

        return D

    def mk_d_spring(self):
        D = np.zeros([4,4],dtype=np.float64)
        D[0,0],D[2,0] =  self.kv, -self.kv
        D[0,2],D[2,2] = -self.kv,  self.kv
        D[1,1],D[3,1] =  self.kh, -self.kh
        D[1,3],D[3,3] = -self.kh,  self.kh

        return D

    # ---------------------------------------------------------
    def mk_visco(self,dof):
        mu = 0.001 # [Pa s]

        if dof == 1:
            D = np.zeros([2,2],dtype=np.float64)
            D[0,0] = mu
            D[1,1] = mu

        elif dof == 2:
            D = np.zeros([3,3],dtype=np.float64)
            D[2,2] = mu

        elif dof == 3:
            D = np.zeros([5,5],dtype=np.float64)
            D[2,2] = mu
            D[3,3] = mu
            D[4,4] = mu

        return D

    # ---------------------------------------------------------
    def mk_imp(self,dof):
        vs = np.sqrt(self.rmu/self.rho)
        vp = np.sqrt((self.rlambda +2*self.rmu)/self.rho)

        if dof == 1:
            imp = np.diag([self.rho*vs])
        elif dof == 2:
            imp = np.diag([self.rho*vp,self.rho*vs])
        elif dof == 3:
            imp = np.diag([self.rho*vp,self.rho*vs,self.rho*vs])

        return imp
