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
