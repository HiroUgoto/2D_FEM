import numpy as np
import copy,math

seed = 1
np.random.seed(seed)
# random.seed(seed)
# t.manual_seed(seed)
# t.cuda.manual_seed(seed)

minlog10x,maxlog10x = -4,5

# GHEモデル
class GHE:
    def __init__(self,G0,gr,hmax,hbeta=1,c1_0=1,c1_inf=1,c2_0=1,c2_inf=1,alpha=0,beta=0):
        self.G0 = G0
        self.gr = gr
        self.hmax = hmax
        self.hbeta = hbeta
        self.c1_0,self.c1_inf,self.c2_0,self.c2_inf = c1_0,c1_inf,c2_0,c2_inf
        self.alpha,self.beta = alpha,beta

        self.gr_h = gr
        self.G0_h = self.G0

        self.tau,self.gamma,self.dg = 0.0,0.0,0.0

        # reversal point
        self.tau0 = 0.0
        self.gamma0 = 0.0

        # p: positive, n: negative
        self.tau0_p_list = []
        self.gamma0_p_list = []

        self.tau0_n_list = []
        self.gamma0_n_list = []

        # yield stress
        self.tau_y = 0.0
        self.gamma_y = 0.0

        self.lmbda0 = 2.0

        self.skeleton = True


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

    def hysteresis_curve(self,gamma):
        tau = self.tau0 + 2*self.skeleton_curve(0.5*(gamma-self.gamma0),self.G0,self.gr)
        return tau

    def get_h(self,gamma,tau):
        return self.hmax*np.abs(1 - (tau/gamma)/self.G0)**self.hbeta

    def update_yield_stress(self,gamma,tau):
        self.tau_y = np.abs(tau)
        self.gamma_y = np.abs(gamma)

    # -------------------------------------------------------------------------------------- #
    def check_reversal(self,gamma,g0,tau0,dg):
        if dg*(gamma-g0) < 0.0:
            self.tau0 = tau0
            self.gamma0 = g0
            self.update_reversal_points(dg)
            self.skeleton = False

    def update_reversal_points(self,dg):
        tau0_list = [self.tau0]
        gamma0_list = [self.gamma0]

        if dg > 0.0:
            for i,tau0 in enumerate(self.tau0_p_list):
                if tau0 > self.tau0:
                    tau0_list = [self.tau0] + self.tau0_p_list[i:]
                    gamma0_list = [self.gamma0] + self.gamma0_p_list[i:]
                    break
            self.tau0_p_list = tau0_list
            self.gamma0_p_list = gamma0_list
            # print("p",self.tau0_p_list,self.tau_y)

        elif dg < 0.0:
            for i,tau0 in enumerate(self.tau0_n_list):
                if tau0 < self.tau0:
                    tau0_list = [self.tau0] + self.tau0_n_list[i:]
                    gamma0_list = [self.gamma0] + self.gamma0_n_list[i:]
                    break
            self.tau0_n_list = tau0_list
            self.gamma0_n_list = gamma0_list
            # print("n",self.tau0_n_list,self.tau_y)

    def find_hysteresis_curve(self,gamma,dg):
        tau = self.hysteresis_curve(gamma)
        if np.abs(tau) >= self.tau_y or self.skeleton:
            self.skeleton = True
            tau = self.skeleton_curve(gamma,self.G0,self.gr)
            self.update_yield_stress(gamma,tau)
            self.tau0_p_list,self.gamma0_p_list = [],[]
            self.tau0_n_list,self.gamma0_n_list = [],[]
            return tau

        elif len(self.tau0_p_list) <= 0:
            return tau
        elif len(self.tau0_n_list) <= 0:
            return tau

        elif dg > 0.0:
            if tau > self.tau0_p_list[0]:
                self.tau0_p_list.pop(0)
                self.gamma0_p_list.pop(0)
                self.tau0_n_list.pop(0)
                self.gamma0_n_list.pop(0)

                if len(self.tau0_n_list)==0:
                    self.tau0 = -self.tau_y
                    self.gamma0 = -self.gamma_y
                else:
                    self.tau0 = self.tau0_n_list[0]
                    self.gamma0 = self.gamma0_n_list[0]
                # print("pop positive",self.tau0,self.gamma0)

                tau = self.hysteresis_curve(gamma)
                return tau
        elif dg < 0.0:
            if tau < self.tau0_n_list[0]:
                self.tau0_p_list.pop(0)
                self.gamma0_p_list.pop(0)
                self.tau0_n_list.pop(0)
                self.gamma0_n_list.pop(0)

                if len(self.tau0_p_list)==0:
                    self.tau0 = self.tau_y
                    self.gamma0 = self.gamma_y
                else:
                    self.tau0 = self.tau0_p_list[0]
                    self.gamma0 = self.gamma0_p_list[0]
                # print("pop negative",self.tau0,self.gamma0)

                tau = self.hysteresis_curve(gamma)
                return tau
        return tau


    # -------------------------------------------------------------------------------------- #
    def shear(self,strain):
        self.check_reversal(strain,self.gamma,self.tau,self.dg)
        if np.abs(strain-self.gamma)>0.0:
            self.dg = strain-self.gamma
        self.tau = self.find_hysteresis_curve(strain,self.dg)
        self.gamma = strain
        return self.tau

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


class GHE_S(GHE):
    def __init__(self,G0,gr,hmax,hbeta=1,c1_0=1,c1_inf=1,c2_0=1,c2_inf=1,alpha=0,beta=0,gr0bygr=5,GmbyGM=0.25):
        super().__init__(G0,gr,hmax,hbeta,c1_0,c1_inf,c2_0,c2_inf,alpha,beta)
        self.gr0 = gr0bygr * gr
        self.GmbyGM = GmbyGM
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

            lmbda = (2-lmbda0)*(gamma/gamma0)**2 + lmbda0
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
            lmbda = (2-self.lmbda0)/gy**2 * dg*(dg+2*gy) + 2
        else:
            lmbda = 2.0
        tau = self.tau0 + lmbda*self.skeleton_curve(dg/lmbda,self.G0_h,self.gr_h)
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
