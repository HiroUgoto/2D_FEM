# from tkinter import E
import numpy as np
import copy, math, sys

import warnings,traceback,sys
warnings.filterwarnings('error')

seed = 1
np.random.seed(seed)
# random.seed(seed)
# t.manual_seed(seed)
# t.cuda.manual_seed(seed)

minlog10x,maxlog10x = -4,5

S_PARAM = {'c1_0':1.0,'c1_inf':0.17,'c2_0':0.8,'c2_inf':2.5,'alpha':2.86,'beta':3.229}
# S_PARAM = {'c1_0':1.0,'c1_inf':0.5,'c2_0':,'c2_inf':2.5,'alpha':2.86,'beta':3.229}


# GHEモデル
class GHE:
    def __init__(self,G0,gr,hmax=0.2,id=None,hbeta=1,c1_0=1,c1_inf=1,c2_0=1,c2_inf=1,alpha=0,beta=0):
        self.G0 = G0
        self.gr = gr
        self.hmax = hmax
        self.hbeta = hbeta
        # self.e0 = e0
        self.c1_0,self.c1_inf,self.c2_0,self.c2_inf = c1_0,c1_inf,c2_0,c2_inf
        self.alpha,self.beta = alpha,beta

        self.gr_h = gr
        self.G0_h = self.G0

        self.tau,self.gamma,self.dg = 0.0,0.0,0.0

        if id is not None:
            self.id = id
        else:
            self.id = -1

        # reversal point
        self.tau0 = 0.0
        self.gamma0 = 0.0

        # p: positive, n: negative
        self.tau0_p_list = [0.0]
        self.gamma0_p_list = [0.0]

        self.tau0_n_list = [0.0]
        self.gamma0_n_list = [0.0]

        # yield stress
        self.tau_y = 0.0
        self.gamma_y = 0.0

        self.lmbda0 = 2.0

        self.skeleton = True

        # h params
        self.xtop = 10
        gamma_top = self.xtop*self.gr
        self.htop = self.get_h(gamma_top,self.skeleton_curve(gamma_top,self.G0,self.gr))
        self.h_rate = 0.05
        self.hmin = 0.5*self.hmax

        # identity
        self.Z3 = np.zeros([3,3])
        self.I3 = np.eye(3)
        self.Dijkl = np.einsum('ij,kl->ijkl',self.I3,self.I3)
        self.Dikjl = np.einsum('ij,kl->ikjl',self.I3,self.I3)
        self.Diljk = np.einsum('ij,kl->ilkj',self.I3,self.I3)
        self.epsilon = sys.float_info.epsilon


    def re_init(self):
        self.G0_h,self.gr_h = self.G0,self.gr
        self.tau,self.gamma,self.dg = 0.0,0.0,0.0
        self.tau0,self.gamma0 = 0,0
        self.tau0_p_list = []
        self.gamma0_p_list = []
        self.tau0_n_list = []
        self.gamma0_n_list = []
        self.tau_y,self.gamma_y = 0,0
        self.lmbda0 = 2.0
        self.skeleton = True

    class StateParameters: #FEM用追加(stress_shiftについて確認する,stress_shiftをカットしている．)
        # def __init__(self,strain,stress,dstrain,dstress,stress_shift,ef1=False,ef2=False):
        def __init__(self,strain,stress,dstrain,dstress,ef1=False,ef2=False):
            self.strain = np.copy(strain)
            self.stress = np.copy(stress)
            self.dstress = np.copy(dstress)
            self.dstrain = np.copy(dstrain)
            self.pmin = 1.0

            self.set_stress_variable()
            self.set_stress_increment()

            # self.elastic_flag1 = ef1
            # self.elastic_flag2 = ef2

        def vector_to_matrix(self,vec):
            mat = np.array([[vec[0],vec[3],vec[5]],
                        [vec[3],vec[1],vec[4]],
                        [vec[5],vec[4],vec[2]]])
            return mat

        def set_stress_variable(self):
            self.p = np.trace(self.stress)/3
            self.sij = self.stress - self.p*np.eye(3)
            if self.p == 0:
                self.rij = 0
            else:
                self.rij = self.sij / self.p
            self.R = np.sqrt(1.5*np.square(self.rij).sum())

        def set_stress_increment(self):
            stress = self.stress + self.dstress
            p = np.trace(stress)/3
            self.dp = p - self.p

    # -------------------------------------------------------------------------------------- #
    def skeleton_curve(self,gamma,G0,gr):
        x = np.abs(gamma/gr)
        if x > 0.0:
            c1 = (self.c1_0+self.c1_inf + (self.c1_0-self.c1_inf)*np.cos(np.pi/(self.alpha/x +1)))/2
            c2 = (self.c2_0+self.c2_inf + (self.c2_0-self.c2_inf)*np.cos(np.pi/(self.beta/x +1)))/2
        else:
            c1,c2 = self.c1_0,self.c2_0
        tau = G0 * gamma /(1/c1 + x/c2)
        return tau

    def hysteresis_curve(self,gamma,gamma0=None,tau0=None):
        if gamma0 is None:
            gamma0 = self.gamma0
        if tau0 is None:
            tau0 = self.tau0
        tau = tau0 + 2*self.skeleton_curve(0.5*(gamma-gamma0),self.G0,self.gr)
        return tau

    def get_h(self, gamma, tau):
        def h1(G):
            return self.hmax*np.max(1-G/self.G0,0)**self.hbeta
        def h2(x):
            return -self.h_rate*(np.log10(x)-1)+self.htop

        x = gamma/self.gr
        if type(x) is np.ndarray:
            h = np.zeros_like(x)
            G = tau/gamma
            h[x<=self.xtop] = h1(G[x<=self.xtop])
            h[x>self.xtop] = np.maximum(h2(x[x>self.xtop]),self.hmin)
        else:
            if x <= self.xtop:
                if abs(gamma) > 0.0:
                    h = h1(tau/gamma)
                else:
                    h = 0.0
            else:
                h = max(h2(x),self.hmin)
        return h


    def update_yield_stress(self,gamma,tau):
        self.tau_y = np.abs(tau)
        self.gamma_y = np.abs(gamma)

    # -------------------------------------------------------------------------------------- #
    def check_reversal(self,gamma,g0,tau0,dg):
        if dg*(gamma-g0) < -self.epsilon:
            self.tau0 = tau0
            self.gamma0 = g0
            self.update_reversal_points(gamma-g0)
            self.skeleton = False
            # print("reverse",gamma-g0,dg,gamma,tau0,self.tau0_p_list,self.tau0_n_list)

    def update_reversal_points(self,dg):
        gamma0_list = [self.gamma0]
        tau0_list = [self.tau0]

        if dg > self.epsilon:
            for i,tau0 in enumerate(self.tau0_n_list):
                if tau0 < self.tau0:
                    tau0_list = [self.tau0] + self.tau0_n_list[i:]
                    gamma0_list = [self.gamma0] + self.gamma0_n_list[i:]
                    break
            self.tau0_n_list = tau0_list
            self.gamma0_n_list = gamma0_list
            # print("n",self.tau0_n_list,self.tau_y)

        elif dg < -self.epsilon:
            for i,tau0 in enumerate(self.tau0_p_list):
                if tau0 > self.tau0:
                    tau0_list = [self.tau0] + self.tau0_p_list[i:]
                    gamma0_list = [self.gamma0] + self.gamma0_p_list[i:]
                    break
            self.tau0_p_list = tau0_list
            self.gamma0_p_list = gamma0_list
            # print("p",self.tau0_p_list,self.tau_y)

        # if len(self.tau0_n_list) < len(self.tau0_p_list):
        #     n = len(self.tau0_n_list)
        #     for i in range(n,len(self.tau0_p_list)):
        #         self.tau0_n_list += [-self.tau_y] 
        #         self.gamma0_n_list += [-self.gamma_y] 

        # elif len(self.tau0_n_list) > len(self.tau0_p_list):
        #     n = len(self.tau0_p_list)
        #     for i in range(n,len(self.tau0_n_list)):
        #         self.tau0_p_list += [self.tau_y] 
        #         self.gamma0_p_list += [self.gamma_y] 


    def find_hysteresis_curve(self,gamma,dg):
        tau = self.find_hysteresis_curve_(gamma,dg)
        self.update_hysteresis_curve(gamma,tau,dg)
        return tau

    def find_hysteresis_curve_(self,gamma,dg):
        if np.abs(gamma) >= self.gamma_y:
            tau = self.skeleton_curve(gamma,self.G0,self.gr)
            # print("A",self.dg,dg,gamma,tau,self.tau0_p_list,self.tau0_n_list)
            return tau

        elif dg > self.epsilon:
            if dg*self.dg < -self.epsilon:
                tau0 = self.tau
                gamma0 = gamma - dg
            else:
                tau0 = -self.tau_y
                gamma0 = -self.gamma_y
            for i in range(len(self.gamma0_p_list)):
                if gamma > self.gamma0_p_list[i]:
                    tau0 = self.tau0_n_list[i]
                    gamma0 = self.gamma0_n_list[i]
                    break
            tau = self.hysteresis_curve(gamma,gamma0,tau0)
            # print("B",self.dg,dg,gamma,tau,self.tau0_p_list,self.tau0_n_list)
            return tau

        elif dg < -self.epsilon:
            if dg*self.dg < -self.epsilon:
                tau0 = self.tau
                gamma0 = gamma - dg
            else:
                tau0 = self.tau_y
                gamma0 = self.gamma_y
            for i in range(len(self.gamma0_n_list)):
                if gamma < self.gamma0_n_list[i]:
                    tau0 = self.tau0_p_list[i]
                    gamma0 = self.gamma0_p_list[i]
                    break
            tau = self.hysteresis_curve(gamma,gamma0,tau0)
            # print("C",self.dg,dg,gamma,tau,self.tau0_p_list,self.tau0_n_list)
            return tau

        else:
            tau = self.hysteresis_curve(gamma)
            # print("D",self.dg,dg,gamma,tau,self.tau0_p_list,self.tau0_n_list)
            return tau

    def update_hysteresis_curve(self,gamma,tau,dg):
        if dg*self.dg < -self.epsilon:
            self.tau0 = self.tau
            self.gamma0 = gamma - dg
            self.update_reversal_points(dg)
            print("reverse",self.id,self.dg,dg,gamma,tau,self.tau0_p_list,self.tau0_n_list)

        if np.abs(gamma) >= self.gamma_y:
            self.update_yield_stress(gamma,tau)
            self.tau0_p_list,self.gamma0_p_list = [self.tau_y],[self.gamma_y]
            self.tau0_n_list,self.gamma0_n_list = [-self.tau_y],[-self.gamma_y]

        elif dg > self.epsilon:
            n = len(self.gamma0_p_list)
            for i in range(n):
                if gamma > self.gamma0_p_list[0]:
                    self.tau0_p_list.pop(0)
                    self.gamma0_p_list.pop(0)
                    self.tau0_n_list.pop(0)
                    self.gamma0_n_list.pop(0)

            if len(self.gamma0_n_list) == 0:
                self.tau0 = -self.tau_y
                self.gamma0 = -self.gamma_y
            else:
                self.tau0 = self.tau0_n_list[0]
                self.gamma0 = self.gamma0_n_list[0]

        elif dg < -self.epsilon:
            n = len(self.gamma0_n_list)
            for i in range(n):
                if gamma < self.gamma0_n_list[0]:
                    self.tau0_p_list.pop(0)
                    self.gamma0_p_list.pop(0)
                    self.tau0_n_list.pop(0)
                    self.gamma0_n_list.pop(0)

            if len(self.gamma0_p_list) == 0:
                self.tau0 = self.tau_y
                self.gamma0 = self.gamma_y
            else:
                self.tau0 = self.tau0_p_list[0]
                self.gamma0 = self.gamma0_p_list[0]

        # print(self.id,dg,gamma,tau,self.tau0_p_list,self.tau0_n_list)


    # -------------------------------------------------------------------------------------- #
    def shear1(self,gamma):
        dg = gamma - self.gamma
        tau = self.find_hysteresis_curve(gamma,dg)
        if np.abs(dg) > self.epsilon:
            self.dg = dg
        self.gamma = gamma
        self.tau = tau
        return self.tau

    def shear_(self,gamma):
        dg = gamma - self.gamma
        tau = self.find_hysteresis_curve_(gamma,dg)
        # print(self.id,dg,gamma,tau)
        return tau

    def cyclic_shear(self,shear_strain,plot=False):
        nstep = len(shear_strain)
        gamma_list,tau_list = [],[]
        for i in range(nstep):
            self.shear(shear_strain[i])
            gamma_list += [self.gamma]
            tau_list += [self.tau]
            # print(gamma,tau,self.gamma0,self.tau0,self.tau_y)
        return np.array(tau_list)


# 前の論文で用いたもの: h-gamma用のgamma_refをパラメータとして用意する
class Prev_IY(GHE):
    def __init__(self,G0,gr,gr_h,hmax,hbeta=1,c1_0=1,c1_inf=1,c2_0=1,c2_inf=1,alpha=0,beta=0):
        super().__init__(G0,gr,hmax,hbeta,c1_0,c1_inf,c2_0,c2_inf,alpha,beta)
        self.gr_h = gr_h

        x = 10.0 ** np.linspace(-5,3,100)
        self.x_for_skeleton = x
        self.h_for_skeleton = 4/np.pi*(1+1/x)*(1 - np.log(1+x)/x) - 2/np.pi

    def update_yield_stress(self,gamma,tau):
        super().update_yield_stress(gamma,tau)
        if np.abs(gamma) > 0.0:
            self.G0_h = np.abs(tau/gamma) * (1+np.abs(gamma)/self.gr_h)

    def hysteresis_curve(self,gamma):
        def skeleton_for_hyst(gamma):
            return self.G0_h * gamma / (1+np.abs(gamma)/self.gr_h)
        tau = self.tau0 + 2*skeleton_for_hyst(0.5*(gamma-self.gamma0))
        return tau

# 修正GHEモデル: h-gammaを数式でモデル化し、
class mod_GHE(GHE):
    def __init__(self,G0,gr,hmax,hbeta=1,c1_0=1,c1_inf=1,c2_0=1,c2_inf=1,alpha=0,beta=0):
        super().__init__(G0,gr,hmax,hbeta,c1_0,c1_inf,c2_0,c2_inf,alpha,beta)
        self.init_hyst()

    def init_hyst(self):
        def skeleton_curve(gamma,G0,gr):
            x = np.abs(gamma/gr)
            idx = x>0
            c1 = np.ones_like(x)*self.c1_0
            c2 = np.ones_like(x)*self.c2_0
            c1[idx] = (self.c1_0+self.c1_inf + (self.c1_0-self.c1_inf)*np.cos(np.pi/(self.alpha/x[idx] +1)))/2
            c2[idx] = (self.c2_0+self.c2_inf + (self.c2_0-self.c2_inf)*np.cos(np.pi/(self.beta/x[idx] +1)))/2
            tau = G0 * gamma /(1/c1 + x/c2)
            return tau

        def hysteresis(gamma,gamma0,tau0,G0_h,gr_h):
            tau = tau0 + 2*skeleton_curve((gamma-gamma0)/2,G0_h,gr_h)
            return tau

        def calculate_h(gamma,tau,gamma0,tau0):
            dW = -np.sum((tau[:-1]+tau[1:])*(gamma[1:]-gamma[:-1])/2, axis=0)
            W = gamma0[0]*tau0[0]/2
            h = 1/(4*np.pi) * dW/W
            return h

        def get_G0_h(gamma0,tau0,gr_h):
            x = np.abs(gamma0/gr_h)
            idx = x>0
            c1 = np.ones_like(x)*self.c1_0
            c2 = np.ones_like(x)*self.c2_0
            c1[idx] = (self.c1_0+self.c1_inf + (self.c1_0-self.c1_inf)*np.cos(np.pi/(self.alpha/x[idx] +1)))/2
            c2[idx] = (self.c2_0+self.c2_inf + (self.c2_0-self.c2_inf)*np.cos(np.pi/(self.beta/x[idx] +1)))/2
            G0_h = np.abs(tau0/gamma0) * (1/c1 + x/c2)
            return G0_h


        n0,n1,n2 = 150,100,1001
        maxlog = np.log10(0.5/self.gr)
        maxlog = max(maxlog10x,maxlog)
        x = 10.0 ** np.linspace(minlog10x,maxlog,n1)[None,:,None]
        gamma0 = x * self.gr

        gamma = np.sin(np.linspace(0,3*np.pi,n0))[:,None,None]
        gamma = gamma*gamma0
        gamma0 = gamma0 * np.ones([n0,n1])[:,:,None]
        gamma0[1:] *= np.sign(np.sign(gamma[1:]-gamma[:-1])+0.5)
        tau0 = skeleton_curve(gamma0,self.G0,self.gr)

        logmin = np.log10(self.gr)
        logmax = logmin + 5
        gr_h = 10.0 ** np.linspace(logmin,logmax,n2)[None,None,:]
        G0_h = get_G0_h(gamma0,tau0,gr_h)
        tau = hysteresis(gamma,gamma0,tau0,G0_h,gr_h)
        k = int(n0/3)
        h_cal = calculate_h(gamma[k:],tau[k:],gamma0,tau0)  # shape: n1,n2
        h = self.get_h(gamma0[0],tau0[0])

        err = np.abs(h_cal - h)
        idx = np.argmin(err,axis=1)
        self.gamma0_for_skeleton = gamma0[0].flatten()
        self.gr_h_for_skeleton = gr_h.flatten()[idx]
        self.G0_h_for_skeleton = G0_h.flatten()[idx]

    def hysteresis_curve(self,gamma):
        tau = self.tau0 + 2*self.skeleton_curve((gamma-self.gamma0)/2,self.G0_h,self.gr_h)
        return tau

    def update_yield_stress(self,gamma,tau):
        def set_hyst_params(gamma,tau):
            def set_G0_h(gamma,tau):
                x = np.abs(gamma/self.gr_h)
                if x > 0.0:
                    c1 = (self.c1_0+self.c1_inf + (self.c1_0-self.c1_inf)*np.cos(np.pi/(self.alpha/x +1)))/2
                    c2 = (self.c2_0+self.c2_inf + (self.c2_0-self.c2_inf)*np.cos(np.pi/(self.beta/x +1)))/2
                    self.G0_h = np.abs(tau/gamma) * (1/c1 + x/c2)
                else:
                    self.G0_h = self.G0
            gamma = np.abs(gamma)
            gamma0 = self.gamma0_for_skeleton
            # G0_h = self.G0_h_for_skeleton
            gr_h = self.gr_h_for_skeleton
            idx1 = np.argmin(np.abs(gamma-gamma0))
            if gamma0[idx1] > gamma:
                idx2 = idx1 - 1
            else:
                idx2 = idx1 + 1
            rate = (gamma-gamma0[idx1])/(gamma0[idx2]-gamma0[idx1])
            # self.G0_h = G0_h[idx1] + rate*(G0_h[idx2]-G0_h[idx1])
            self.gr_h = gr_h[idx1] + rate*(gr_h[idx2]-gr_h[idx1])
            set_G0_h(gamma,tau)

        set_hyst_params(gamma,tau)
        super().update_yield_stress(gamma,tau)


# 修正GHEモデル： 履歴曲線用の骨格曲線をHDモデルにしたパターン
class mod_GHE2(GHE):
    def hysteresis_curve(self,gamma):
        def skeleton_curve(gamma):
            x = np.abs(gamma/self.gr_h)
            return self.G0_h*gamma / (1+x)
        tau = self.tau0 + 2*skeleton_curve((gamma-self.gamma0)/2)
        return tau

    def update_yield_stress(self,gamma,tau):
        def set_hyst_params(gamma,tau):
            def set_gr_h(gamma,tau):
                h = self.get_h(gamma,tau)
                if h > 0.0:
                    logx = np.arctanh(np.pi*h-1) + 0.64
                    x = 10.0 ** logx
                    self.gr_h = gamma/x

            def set_G0_h(gamma,tau):
                x = np.abs(gamma/self.gr_h)
                self.G0_h = np.abs(tau/gamma)*(1+x)
            if gamma > 0.0:
                set_gr_h(gamma,tau)
                set_G0_h(gamma,tau)

        set_hyst_params(np.abs(gamma),np.abs(tau))
        super().update_yield_stress(gamma,tau)


class GHE_S(GHE):
    coef = (1,0)

    def __init__(self,G0,gr,hmax,id=None,hbeta=1,c1_0=1,c1_inf=1,c2_0=1,c2_inf=1,alpha=0,beta=0,gr0bygr=5,GmbyGM=0.1):
        super().__init__(G0,gr,hmax,id,hbeta,c1_0,c1_inf,c2_0,c2_inf,alpha,beta)
        self.gr0 = gr0bygr * gr
        self.GmbyGM = GmbyGM
        self.init_hyst()

        # identity
        self.Z3 = np.zeros([3,3])
        self.I3 = np.eye(3)
        self.Dijkl = np.einsum('ij,kl->ijkl',self.I3,self.I3)
        self.Dikjl = np.einsum('ij,kl->ikjl',self.I3,self.I3)
        self.Diljk = np.einsum('ij,kl->ilkj',self.I3,self.I3)
        # exit()
        # print("aaaa")

    def vector_to_matrix(self,vec):
        mat = np.array([[vec[0],vec[3],vec[5]],
                        [vec[3],vec[1],vec[4]],
                        [vec[5],vec[4],vec[2]]])
        return mat

    def clear_strain(self):
        self.strain = self.Z3.copy()

    def elastic_modulus_E(self,G,nu): #FEM追加(solve_strainでひずみを求めるためのEを計算している．)
        rlambda_coeff = 2*nu/(2-nu)
        rlambda = G*rlambda_coeff
        Ee = rlambda*self.Dijkl + G*(self.Dikjl+self.Diljk)
        return Ee

    def solve_strain(self,stress_mat,E): #FEM追加(isotropic_compressionGでひずみを求めるために使用．σ=Eεからεを求めている．)
        b = stress_mat.flatten()
        A = np.reshape(E,(9,9))
        Ainv = np.linalg.pinv(A)
        x = Ainv @ b
        strain_mat = np.reshape(x,(3,3))
        return strain_mat

    def isotropic_compressionG(self,compression_stress,nstep=10): #FEM追加(自重解析のところ，最後にひずみをクリア),nstep=1000
        dcp = compression_stress / nstep
        dstress_vec = np.array([dcp,dcp,dcp,0,0,0])
        dstress = self.vector_to_matrix(dstress_vec)

        self.stress = np.zeros((3,3))
        self.strain = np.zeros((3,3))

        for i in range(nstep):
            E = self.elastic_modulus_E(G=3.903e7,nu=0.40)
            dstrain = self.solve_strain(dstress,E)

            self.stress += dstress
            self.strain += dstrain

        self.clear_strain()

    def set_stress_variable(self,stress): #FEM追加(応力の定義)
        p = np.trace(stress)/3
        r_stress = (stress-p*self.I3) / p
        R = math.sqrt(1.5)*np.linalg.norm(r_stress)
        return p,R

    def set_strain_variable(self,strain):
        ev = np.trace(strain)
        dev_strain = strain - ev/3.0 * self.I3
        gamma = math.sqrt(2.0/3.0)*np.linalg.norm(dev_strain)
        return ev,gamma

        
    def init_hyst(self):
        def get_G0_h(gamma0):
            G0_hbyGmax = (1-self.GmbyGM)/(1+gamma0/self.gr0) + self.GmbyGM
            return G0_hbyGmax * self.G0

        def set_gr_h(gamma0):
            n = len(gamma0)
            gr_h = []
            log_gr = np.log10(self.gr)
            for i in range(n):
                g0 = gamma0[i]
                G0_h = get_G0_h(g0)
                n2 = 1000
                gr = 10.0 ** np.linspace(log_gr,log_gr+3,n2)
                tau0 = skeleton_curve(g0,self.G0,self.gr)
                tau = skeleton_curve(g0,G0_h,gr)
                err = np.abs(tau-tau0)
                idx = np.argmin(err)
                gr_h += [gr[idx]]
            self.gr_h_for_skeleton = np.array(gr_h)

        def skeleton_curve(gamma,G0,gr):
            x = np.abs(gamma/gr)
            if type(x) is np.ndarray:
                idx = x>0
                c1 = np.ones_like(x)*self.c1_0
                c2 = np.ones_like(x)*self.c2_0
                c1[idx] = (self.c1_0+self.c1_inf + (self.c1_0-self.c1_inf)*np.cos(np.pi/(self.alpha/x[idx] +1)))/2
                c2[idx] = (self.c2_0+self.c2_inf + (self.c2_0-self.c2_inf)*np.cos(np.pi/(self.beta/x[idx] +1)))/2
            else:
                if x > 0.0:
                    c1 = (self.c1_0+self.c1_inf + (self.c1_0-self.c1_inf)*np.cos(np.pi/(self.alpha/x +1)))/2
                    c2 = (self.c2_0+self.c2_inf + (self.c2_0-self.c2_inf)*np.cos(np.pi/(self.beta/x +1)))/2
                else:
                    c1,c2 = self.c1_0,self.c2_0
            tau = G0 * gamma /(1/c1 + x/c2)
            return tau

        def hysteresis(gamma,gamma0,tau0,lmbda0):
            gr_h = self.gr_h_for_skeleton[None,:,None]
            G0_h = get_G0_h(np.abs(gamma0))
            x = gamma/gamma0

            # fx = (self.coef[0]*x**2 + self.coef[1]*x**4)/sum(self.coef)
            fx = sum([c*x**(2*i+2) for i,c in enumerate(self.coef)])/sum(self.coef)
            lmbda = (2-lmbda0)*fx + lmbda0
            tau = tau0 + lmbda*skeleton_curve((gamma-gamma0)/lmbda,G0_h,gr_h)
            return tau

        def calculate_h(gamma,tau,gamma0,tau0):
            dW = -np.sum((tau[:-1]+tau[1:])*(gamma[1:]-gamma[:-1])/2, axis=0)
            W = gamma0[0]*tau0[0]/2
            h = 1/(4*np.pi) * dW/W
            return h


        n0,n1,n2 = 150,100,501
        maxlog = np.log10(0.5/self.gr)
        maxlog = max(maxlog10x,maxlog)
        x = 10.0 ** np.linspace(minlog10x,maxlog,n1)[None,:,None]
        gamma0 = x * self.gr

        set_gr_h(gamma0.flatten())

        lmbda0 = np.linspace(1,2,n2)[None,None,:]
        gamma = np.sin(np.linspace(0,3*np.pi,n0))[:,None,None]
        gamma = gamma*gamma0
        gamma0 = gamma0 * np.ones([n0,n1])[:,:,None]
        gamma0[1:] *= np.sign(np.sign(gamma[1:]-gamma[:-1])+0.5)
        tau0 = skeleton_curve(gamma0,self.G0,self.gr)

        tau = hysteresis(gamma,gamma0,tau0,lmbda0)
        k = int(n0/3)
        h_cal = calculate_h(gamma[k:],tau[k:],gamma0,tau0)  # shape: n1,n2
        h = self.get_h(gamma0[0],tau0[0])

        err = np.abs(h_cal - h)
        idx = np.argmin(err,axis=1)
        self.gamma0_for_skeleton = gamma0[0].flatten()
        self.lmbda0_for_skeleton = lmbda0.flatten()[idx]

    def hysteresis_curve(self,gamma,gamma0=None,tau0=None):
        if gamma0 is None:
            gamma0 = self.gamma0
        if tau0 is None:
            tau0 = self.tau0

        dg = gamma - gamma0
        if np.abs(dg)>0.0:
            gy = self.gamma_y * np.sign(-dg)
        else:
            gy = self.gamma_y
        if np.abs(self.gamma_y) > 0:
            x = dg/gy + 1
            # fx = (self.coef[0]*x**2 + self.coef[1]*x**4)/sum(self.coef)
            fx = sum([c*x**(2*i+2) for i,c in enumerate(self.coef)])/sum(self.coef)
            lmbda0 = max(self.lmbda0,1)
            lmbda = (2-lmbda0)*fx + lmbda0
            # lmbda = (2-self.lmbda0)*(self.coef[0]*x**2 + self.coef[1]*x**4) + self.lmbda0
            # lmbda = (2-self.lmbda0)*x**self.degree + self.lmbda0
            # lmbda = (2-self.lmbda0)/gy**self.degree * dg*(dg+2*gy) + 2
        else:
            lmbda = 2.0
        tau = tau0 + lmbda*self.skeleton_curve(dg/lmbda,self.G0_h,self.gr_h)
        return tau

    def update_yield_stress(self,gamma,tau):
        def set_hyst_params(gamma):
            gamma = np.abs(gamma)
            gamma0 = self.gamma0_for_skeleton
            lmbda0 = self.lmbda0_for_skeleton
            gr_h = self.gr_h_for_skeleton
            G0_hbyGmax = (1-self.GmbyGM)/(1+gamma0/self.gr0) + self.GmbyGM
            G0_h = G0_hbyGmax * self.G0

            if gamma < gamma0[0]:
                self.lmbda0 = 2.0
                self.G0_h = self.G0
                self.gr_h = self.gr
            else:
                idx1 = np.argmin(np.abs(np.abs(gamma-gamma0)))
                if gamma0[idx1] > gamma:
                    idx2 = idx1 - 1
                else:
                    idx2 = idx1 + 1
                rate = (gamma-gamma0[idx1])/(gamma0[idx2]-gamma0[idx1])
                self.lmbda0 = lmbda0[idx1] + rate * (lmbda0[idx2]-lmbda0[idx1])
                self.G0_h = G0_h[idx1] + rate * (G0_h[idx2]-G0_h[idx1])
                self.gr_h = gr_h[idx1] + rate * (gr_h[idx2]-gr_h[idx1])
        set_hyst_params(gamma)
        super().update_yield_stress(gamma,tau)


class myGHE_S(GHE_S):
    # coef = (0,0,0,1)
    coef = (1,)
    # coef = (1,1)


class special_GHE_S(GHE):
    coef = (1,)

    # def __init__(self,G0,gr,hmax,hbeta=1,c1_0=1,c1_inf=1,c2_0=1,c2_inf=1,alpha=0,beta=0,gr0bygr=5,GmbyGM=0.18):
    def __init__(self,G0,gr,hmax,hbeta=1,c1_0=1,c1_inf=1,c2_0=1,c2_inf=1,alpha=0,beta=0,gr0bygr=5,GmbyGM=0.1):
        super().__init__(G0,gr,hmax,hbeta,c1_0,c1_inf,c2_0,c2_inf,alpha,beta)
        self.gr0 = gr0bygr * gr
        self.GmbyGM = GmbyGM
        self.xmax = 0
        self.rh = 500
        self.init_hyst()

    def init_hyst(self):
        def get_G0_h(gamma0):
            G0_hbyGmax = (1-self.GmbyGM)/(1+gamma0/self.gr0) + self.GmbyGM
            return G0_hbyGmax * self.G0

        def set_gr_h(gamma0):
            n = len(gamma0)
            gr_h = []
            log_gr = np.log10(self.gr)
            for i in range(n):
                g0 = gamma0[i]
                G0_h = get_G0_h(g0)
                n2 = 1000
                gr = 10.0 ** np.linspace(log_gr,log_gr+3,n2)
                tau0 = skeleton_curve(g0,self.G0,self.gr)
                tau = skeleton_curve(g0,G0_h,gr)
                err = np.abs(tau-tau0)
                idx = np.argmin(err)
                gr_h += [gr[idx]]
            self.gr_h_for_skeleton = np.array(gr_h)

        def skeleton_curve(gamma,G0,gr):
            x = np.abs(gamma/gr)
            if type(x) is np.ndarray:
                idx = x>0
                c1 = np.ones_like(x)*self.c1_0
                c2 = np.ones_like(x)*self.c2_0
                c1[idx] = (self.c1_0+self.c1_inf + (self.c1_0-self.c1_inf)*np.cos(np.pi/(self.alpha/x[idx] +1)))/2
                c2[idx] = (self.c2_0+self.c2_inf + (self.c2_0-self.c2_inf)*np.cos(np.pi/(self.beta/x[idx] +1)))/2
            else:
                if x > 0.0:
                    c1 = (self.c1_0+self.c1_inf + (self.c1_0-self.c1_inf)*np.cos(np.pi/(self.alpha/x +1)))/2
                    c2 = (self.c2_0+self.c2_inf + (self.c2_0-self.c2_inf)*np.cos(np.pi/(self.beta/x +1)))/2
                else:
                    c1,c2 = self.c1_0,self.c2_0
            tau = G0 * gamma /(1/c1 + x/c2)
            return tau

        def hysteresis(gamma,gamma0,tau0,lmbda0):
            gr_h = self.gr_h_for_skeleton[None,:,None]
            G0_h = get_G0_h(np.abs(gamma0))
            x = gamma/gamma0

            # fx = (self.coef[0]*x**2 + self.coef[1]*x**4)/sum(self.coef)
            fx = np.zeros_like(x)
            fx[x<0] = sum([c*x[x<0]**(2*i+2) for i,c in enumerate(self.coef)])/sum(self.coef)
            fx[x>=0] = (self.rh/(self.xmax+self.rh) * x[x>=0])**2
            lmbda = (2-lmbda0)*fx + lmbda0
            tau = tau0 + lmbda*skeleton_curve((gamma-gamma0)/lmbda,G0_h,gr_h)
            return tau

        def calculate_h(gamma,tau,gamma0,tau0):
            dW = -np.sum((tau[:-1]+tau[1:])*(gamma[1:]-gamma[:-1])/2, axis=0)
            W = gamma0[0]*tau0[0]/2
            h = 1/(4*np.pi) * dW/W
            return h


        n0,n1,n2 = 150,100,501
        maxlog = np.log10(0.5/self.gr)
        maxlog = max(maxlog10x,maxlog)
        x = 10.0 ** np.linspace(minlog10x,maxlog,n1)[None,:,None]
        gamma0 = x * self.gr

        set_gr_h(gamma0.flatten())

        lmbda0 = np.linspace(1,2,n2)[None,None,:]
        gamma = np.sin(np.linspace(0,3*np.pi,n0))[:,None,None]
        gamma = gamma*gamma0
        gamma0 = gamma0 * np.ones([n0,n1])[:,:,None]
        gamma0[1:] *= np.sign(np.sign(gamma[1:]-gamma[:-1])+0.5)
        tau0 = skeleton_curve(gamma0,self.G0,self.gr)

        tau = hysteresis(gamma,gamma0,tau0,lmbda0)
        k = int(n0/3)
        h_cal = calculate_h(gamma[k:],tau[k:],gamma0,tau0)  # shape: n1,n2
        h = self.get_h(gamma0[0],tau0[0])

        err = np.abs(h_cal - h)
        idx = np.argmin(err,axis=1)
        self.gamma0_for_skeleton = gamma0[0].flatten()
        self.lmbda0_for_skeleton = lmbda0.flatten()[idx]

    def hysteresis_curve(self,gamma):
        dg = gamma - self.gamma0
        if np.abs(dg)>0.0:
            gy = self.gamma_y * np.sign(-dg)
        else:
            gy = self.gamma_y
        if np.abs(self.gamma_y) > 0:
            x = dg/gy + 1
            # fx = (self.coef[0]*x**2 + self.coef[1]*x**4)/sum(self.coef)
            if x<0:
                fx = sum([c*x**(2*i+2) for i,c in enumerate(self.coef)])/sum(self.coef)
            else:
                fx = (self.rh/(self.xmax+self.rh) * x)**2
            lmbda0 = max(self.lmbda0,1.05)
            lmbda = (2-lmbda0)*fx + lmbda0
            # lmbda = (2-self.lmbda0)*(self.coef[0]*x**2 + self.coef[1]*x**4) + self.lmbda0
            # lmbda = (2-self.lmbda0)*x**self.degree + self.lmbda0
            # lmbda = (2-self.lmbda0)/gy**self.degree * dg*(dg+2*gy) + 2
        else:
            lmbda = 2.0
        tau = self.tau0 + lmbda*self.skeleton_curve(dg/lmbda,self.G0_h,self.gr_h)
        return tau

    def update_yield_stress(self,gamma,tau):
        def get_G0_h(gamma):
            G0_hbyGmax = (1-self.GmbyGM)/(1+gamma/self.gr0) + self.GmbyGM
            return G0_hbyGmax * self.G0

        def set_hyst_params(gamma):
            gamma = np.abs(gamma)
            gamma0 = self.gamma0_for_skeleton
            lmbda0 = self.lmbda0_for_skeleton
            gr_h = self.gr_h_for_skeleton
            G0_hbyGmax = (1-self.GmbyGM)/(1+gamma0/self.gr0) + self.GmbyGM
            G0_h = G0_hbyGmax * self.G0
            if gamma < gamma0[0]:
                self.lmbda0 = 2.0
                self.G0_h = self.G0
                self.gr_h = self.gr
            else:
                idx1 = np.argmin(np.abs(np.abs(gamma-gamma0)))
                if gamma0[idx1] > gamma:
                    idx2 = idx1 - 1
                else:
                    idx2 = idx1 + 1
                rate = (gamma-gamma0[idx1])/(gamma0[idx2]-gamma0[idx1])
                self.lmbda0 = lmbda0[idx1] + rate * (lmbda0[idx2]-lmbda0[idx1])
                self.gr_h = gr_h[idx1] + rate * (gr_h[idx2]-gr_h[idx1])
                # self.G0_h = G0_h[idx1] + rate * (G0_h[idx2]-G0_h[idx1])
                self.G0_h = get_G0_h(gamma)
                self.xmax = gamma/self.gr
        set_hyst_params(gamma)
        super().update_yield_stress(gamma,tau)



class Multi:
    def __init__(self,models=[mod_GHE(50e6,1e-3,0.2,**S_PARAM)]):
        self.models = models
        self.nmodel = len(models)

    def shear(self,strain):
        # tau = np.array([m.shear(strain) for m in self.models]).mean()
        tau = sum([m.shear(strain) for m in self.models])/self.nmodel
        return tau

    def cyclic_shear(self,shear_strain):
        tau = np.zeros([self.nmodel,len(shear_strain)])
        for i,m in enumerate(self.models):
            tau[i] = m.cyclic_shear(shear_strain)
        return tau.mean(axis=0)



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import setfig
    setfig.init({'figure.dpi':200,'figure.figsize':setfig.Square*0.8,'mathtext.fontset':'cm'})

    xa_list = [10,100,1000,10000]
    for xa in xa_list:
        ncycle = 3
        div = 100
        cycle = np.linspace(0,ncycle*2*np.pi,ncycle*div)
        gamma = xa * np.sin(cycle)

        args = {'G0':1,'gr':1,'hmax':0.2}
        args.update(S_PARAM)
        model1 = mod_GHE(**args)
        model2 = GHE_S(**args)
        model3 = special_GHE_S(**args)

        tau1 = model1.cyclic_shear(gamma)[div:]
        tau2 = model2.cyclic_shear(gamma)[div:]
        tau3 = model3.cyclic_shear(gamma)[div:]

        gamma = gamma[div:]
        gr = args['gr']

        fig,ax = plt.subplots()
        ax.set_xlabel(r"正規化ひずみ $x$")
        ax.set_ylabel(r'正規化せん断応力 $y$', labelpad=4.0)
        ax.plot(gamma,tau1,label='修正GHE')
        ax.plot(gamma,tau2,label='GHE-S')
        ax.plot(gamma,tau3,label='修正GHE-S')
        # ax.set_yticks([])
        ax.legend()
        fig.savefig(f'ghe_{xa:d}.png')
        plt.close(fig)