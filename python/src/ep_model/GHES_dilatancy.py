import numpy as np
import math,sys
import matplotlib.pyplot as plt
try:
    import GHES_hd_FEM as hd_FEM
except:
    from . import GHES_hd_FEM as hd_FEM
import warnings,traceback
warnings.filterwarnings('error')

def devortic(vector):
    vec = vector.copy() ##そのままひずみ・応力を返している
    vec[:3] -= vec[:3].mean() ##偏差ひずみ・偏差応力を返している
    v2 = vec**2
    v2[3:] *= 2
    return np.sqrt(v2.sum())


class Martin_Finn_Seed:
    def __init__(self,p0,k2,M=1.25,a=0.7,m=0.43,n=0.62,c1=0.8,c2=0.5):
        self.p0 = p0
        self.M = M/np.sqrt(3)
        # self.M = M
        self.phase_boundary = self.M*a
        self.k2 = k2
        self.m = m
        self.n = n
        self.c1 = c1
        self.c2 = c2

        self.init()

    def init(self):
        self.ev = 0.0
        self.temp_dev = 0.0
        self.prevp = self.p0
        self.prevgamma = 0.0
        self.prevtau = 0.0
        self.maxgamma = 0.0
        self.maxp = 0.0
        self.pmin = 100
        self.mode = 0  # 0:通常状態, 1:載荷(変相線外), 2:除荷(変相線外)
        self.m2 = 0.0
        self.cum_du = 0.0
        self.gamma_list = []
        self.p_list = []


    def dev(self, gamma):
        dev = gamma*self.c1*np.exp(-self.c2*self.ev/gamma)
        self.temp_dev = dev
        return dev

    def get_du(self,gamma):
        Er = self.prevp**(1-self.m)/(self.m*self.k2*self.p0**(self.n-self.m))
        du = Er*self.dev(gamma)
        # print(f'Er:{Er:.1e}, dev:{self.temp_dev:.1e}, du:{du:.1e}')
        return du * (gamma-self.prevgamma)

    def _shear(self,gamma,tau):
        def mode0(gamma,tau):
            du = self.get_du(gamma)
            if self.prevp-du > self.pmin:
                return self.prevp - du
            else:
                return self.pmin

        def mode1(gamma,tau):
            self.cum_du += self.get_du(gamma)
            p = np.sqrt(self.m2 + (tau/self.M)**2)
            return p

        def mode2(gamma,tau):
            cum_du = self.cum_du
            if self.maxp - cum_du < self.pmin:
                pmin = self.pmin
            else:
                pmin = self.maxp - cum_du
            p = (self.maxp-pmin)/self.maxgamma * gamma + pmin
            dp = p - self.maxp  # <0
            Er = self.prevp**(1-self.m)/(self.m*self.k2*self.p0**(self.n-self.m))
            self.temp_dev = -dp/Er
            self.maxgamma = gamma
            self.maxp = p
            self.cum_du += dp
            if self.cum_du < 0.0:
                self.cum_du = 0.0
            return p

        def to_mode1(gamma,tau,p):
            self.mode = 1
            if np.isclose(p,self.prevp):
                p_cross = self.prevp
            else:
                r = (tau-self.prevtau)/(p-self.prevp)
                p_cross = (self.prevtau-r*self.prevp)/(self.phase_boundary-r)
            tau_cross = self.phase_boundary*p_cross
            self.m2 = p_cross**2 - (tau_cross/self.M)**2
            self.prevp = mode1(gamma,tau)

        gamma = np.abs(gamma)
        tau = np.abs(tau)
        dgamma = gamma - self.prevgamma
        dtau = tau - self.prevtau
        loading = dgamma>0 and dtau>0

        mode = self.mode
        prevp = self.prevp
        if self.mode == 0:
            if loading:
                p = mode0(gamma,tau)
                if tau/p > self.phase_boundary:
                    to_mode1(gamma,tau,p)
                else:
                    self.prevp = p

        elif self.mode == 1:
            if dtau >= 0:
                self.prevp = mode1(gamma,tau)
            else:
                self.mode = 2
                self.maxgamma = self.prevgamma
                self.maxp = self.prevp
                self.prevp = mode2(gamma,tau)

        elif self.mode == 2:
            if loading:
                self.mode = 0
                p = mode0(gamma,tau)
            else:
                p = mode2(gamma,tau)

            if (tau/p > self.phase_boundary) and (dtau>0):
                to_mode1(gamma,tau,p)
            else:
                self.prevp = p

        dp = self.prevp - prevp
        # print(f'mode:{mode}, eta:{tau/self.prevp}')
        # if (tau/self.prevp > self.phase_boundary) and (self.mode == 0):
        # if (self.mode==2) and (dp>0):
        if (self.mode==1) and (dtau<0):
            print(mode,self.mode)

        self.prevgamma = gamma
        self.prevtau = tau
        self.ev += self.temp_dev
        self.gamma_list.append(gamma)
        self.p_list.append(self.prevp)
        return self.prevp

    def shear(self,strain_vec,stress_vec):
        e = strain_vec.copy()
        e[:3] -= strain_vec[:3].mean()
        gamma = np.sqrt(2 * (e**2).sum())  # gamma for sheer
        s = stress_vec.copy()
        s[:3] -= stress_vec[:3].mean()
        tau = np.sqrt((s**2).sum())
        return self._shear(gamma,tau)


class Martin_Finn_Seed2:
    def __init__(self,p0,k2=1e-4,M=1.5,a=0.4,m=0.43,n=0.62,c1=1.0,c2=0.4):
        self.p0 = p0
        # self.M = M/np.sqrt(3)
        self.M = M
        # self.phase_boundary = self.M*a
        self.a_phase = (1-a)/a**2/self.M**2/self.p0
        self.k2 = k2
        self.m = m
        self.n = n
        self.c1 = c1
        self.c2 = c2
        self.c2 = 0.01
        self.c2 = 0.5

        self.init()

    def init(self):
        self.ev = 0.0
        self.temp_dev = 0.0
        self.p_cycle = self.p0
        self.p_max_cycle = self.p0
        self.prevp = self.p0
        self.prevgamma = 0.0
        self.prevtau = 0.0
        self.taumax = 0.0
        self.pmax = 0.0
        self.pmin = 1e3
        self.gamma_min = 1e-6
        self.mode = 0  # 0:通常状態, 1:載荷(変相線外), 2:除荷(変相線外)
        self.m2 = 0.0
        self.cum_du = 0.0
        self.gamma_list = []
        self.p_list = []
        self.mode_list = []

    def phase_boundary(self,tau):
        return self.a_phase*tau**2 + tau/self.M

    def dev(self, gamma):
        # gamma = gamma**0.5
        # gamma_min = 0.05
        # gamma = np.sqrt(gamma**2+gamma_min**2)
        # dev = gamma*self.c1*np.exp(-self.c2*self.ev/gamma) * self.dgamma
        dev = 0.1*self.c1*np.exp(-self.c2*self.ev/gamma) * self.dgamma  #体積ひずみ増分の式(3.11)
        self.temp_dev = dev
        return dev

    def Er(self): #膨潤係数の式(3.12)
        # Er = self.prevp**(1-self.m)/(self.m*self.k2*self.p0**(self.n-self.m))
        Er = self.p_cycle**(1-self.m)/(self.m*self.k2*self.p0**(self.n-self.m))
        return Er

    def get_du(self,gamma): #間隙水圧の計算式(3.13)
        du = 0.0
        # if gamma > self.gamma_min:
        #     du = self.Er()*self.dev(gamma)
        # else:
        #     du = 0.0
        return du

    def gamma_process(self,gamma):
        if self.dgamma > 0:
            du = self.get_du(gamma)
            self.p_cycle = max(self.p_cycle-du, self.pmin)
            self.ev += self.temp_dev

    def mode0(self,tau):
        self.mode = 0
        return self.p_cycle

    def mode1(self,tau):
        self.mode = 1
        p = np.sqrt(self.m2 + (tau/self.M)**2)
        self.taumax = tau
        self.pmax = p
        self.p_max_cycle = max(p,self.p_max_cycle)
        return p

    def mode2(self,tau):
        if self.mode != 2:
            self.p_max_cycle = self.prevp
        self.mode = 2
        p = (self.pmax-self.p_cycle)/self.taumax * tau + self.p_cycle
        self.taumax = tau
        self.pmax = p
        return p

    def tau_process(self,tau):
        mode = self.mode
        if mode == 0:
            p = self.mode0(tau)
            if p < self.phase_boundary(tau):
                p_cross = self.prevp
                tau_cross = self.prevtau
                self.m2 = p_cross**2 - (tau_cross/self.M)**2
                self.prevp = self.mode1(tau)
            else:
                self.prevp = p

        elif mode == 1:
            if self.dtau >= 0:
                self.prevp = self.mode1(tau)
            else:
                self.prevp = self.mode2(tau)

        elif mode == 2:
            if self.dtau < 0:
                self.prevp = self.mode2(tau)
            else:
                if self.prevp < self.phase_boundary(self.prevtau):
                    self.m2 = self.prevp**2 - (self.prevtau/self.M)**2
                    self.prevp = self.mode1(tau)
                else:
                    p = self.mode0(tau)
                    if p < self.phase_boundary(tau):
                        p_cross = self.prevp
                        tau_cross = self.prevtau
                        self.m2 = p_cross**2 - (tau_cross/self.M)**2
                        self.prevp = self.mode1(tau)
                    else:
                        self.prevp = p


    def _shear(self,gamma,tau):
        gamma,tau = np.abs(gamma),np.abs(tau)
        self.dgamma = gamma - self.prevgamma
        self.dtau = tau - self.prevtau
        self.gamma_process(gamma)
        self.tau_process(tau)
        self.prevgamma = gamma
        self.prevtau = tau
        self.gamma_list.append(gamma)
        self.p_list.append(self.prevp)
        self.mode_list.append(self.mode)
        return self.prevp

    def shear(self,strain_vec,stress_vec):
        gamma = devortic(strain_vec)
        tau = devortic(stress_vec)
        return self._shear(gamma,tau)


class Martin_Finn_Seed3(Martin_Finn_Seed2):
    def __init__(self,p0,k2=1e-4,M=1.25,a=0.2,m=0.43,n=0.62,c1=1.0,c2=0.4):
        return super().__init__(p0,k2,M,a,m,n,c1,c2)

    def init(self):
        super().init()
        self.prevdg = 0.0

    def dev(self,gamma):
        dev = gamma*self.c1*np.exp(-self.c2*self.ev/gamma)
        self.temp_dev = dev
        return dev

    def gamma_process(self, gamma):
        if (self.dgamma<0) and (self.prevdg>0):
            du = self.get_du(gamma)
            self.p_cycle = max(self.p_cycle-du,self.pmin)
            self.ev += self.temp_dev
        if abs(self.dgamma)>0:
            self.prevdg = self.dgamma

    def mode0(self,tau):
        self.mode = 0
        if (self.dtau>0) or (np.abs(self.prevtau)<=0.0):
            return self.prevp
        else:
            p = (self.prevp-self.p_cycle)/self.prevtau * tau + self.p_cycle
            return p


if __name__ == '__main__':
    import hd_FEM
    import setfig
    setfig.init({'axes.labelsize':20,})

    p0 = 50e3
    ncycle = 50
    div = 100
    amp = 0.15/2

    ncycle = 5-0.25
    amp = 0.05

    cycle = np.linspace(0,ncycle,math.floor(ncycle*div))
    ea = amp * np.sin(cycle*2*np.pi)
    # ea[:int(len(ea))] *= np.linspace(0,1,int(len(ea)))
    ea[:int(len(ea))] *= np.linspace(0,1,int(len(ea)))
    # ea += np.linspace(0,2*amp,int(len(ea)))

    # fig,ax = plt.subplots()
    # ax.plot(ea) #label
    # plt.show()
    # plt.close(fig)
    gamma = np.sqrt(3)*ea

    phi = 30

    M = 6*np.sin(np.radians(phi))/(3-np.sin(np.radians(phi)))
    # tauf = p0*M/np.sqrt(3) / 2
    tauf = M*p0*np.cos(np.radians(phi))/2
    gr,hmax = 2.13e-4,0.2 #値いじる前
    # gr,hmax = 5e-4,0.2
    G0 = tauf/(2.5*gr)
    # M = 1.5
    # tauf = p0*M/np.sqrt(3) / 2
    # gr,hmax = 1e-3,0.2
    # G0 = tauf/(2.5*gr)
    ghe = hd_FEM.special_GHE_S(G0=G0,gr=gr,hmax=hmax,**hd_FEM.S_PARAM)
    tau = ghe.cyclic_shear(gamma)
    tau_max = tau.max()*1.2

    q = np.sqrt(3)*tau
    # model = Martin_Finn_Seed(p0,k2=1e-7)
    # model = Martin_Finn_Seed2(p0,k2=1e-4)
    # model = Martin_Finn_Seed2(p0,k2=1e-5)
    model = Martin_Finn_Seed3(p0,M=M,k2=3e-4)
    p = []
    for i in range(len(gamma)):
        p.append(model._shear(gamma[i],tau[i]))
    p = np.array(p)/1000
    p0 /= 1000
    tau /= 1000
    mode = np.array(model.mode_list)

    tau_phase = np.linspace(-tau_max,tau_max,200)
    p_phase = model.phase_boundary(np.abs(tau_phase)) / 1000
    p_phase[p_phase>p0] = np.nan
    p_crit = np.abs(tau_phase)/M / 1000
    tau_phase /= 1000

    fig,ax = plt.subplots(figsize=setfig.Short)
    ax.set_xlabel(f"$p'$")
    ax.set_ylabel(r'$q$', labelpad=4.0)
    lx = np.array([-p0,0,p0])*0.6
    lx1 = lx*0.7
    lx2 = lx*0.6
    # l1 = lx1*model.phase_boundary
    l2 = lx2*model.M
    # lx1 = np.abs(lx1)
    lx2 = np.abs(lx2)
    # ax.plot(p/1000,tau/1000,lw=0.5,marker='o',markersize=1,antialiased=False)
    ax.plot(p_crit,tau_phase,c='red',label='破壊線')
    ax.plot(p_phase,tau_phase,c='blue',label='変相線')
    ax.plot(p,tau,c='black',label='モデル')
    # ax.scatter(p[mode==0],tau[mode==0],s=1,c='blue',label='mode:0')
    # ax.scatter(p[mode==1],tau[mode==1],s=1,c='orange',label='mode:1')
    # ax.scatter(p[mode==2],tau[mode==2],s=1,c='green',label='mode:2')
    # ax.plot(lx1,l1,lw=0.5,c='black')
    # ax.plot(lx2,l2,lw=0.5,c='black')
    # ax.legend()
    ax.set_xlim(0,p.max()*1.2)
    ax.set_ylim(tau.min()*1.3,tau.max()*1.3)
    setfig.legend(ax,'upper right')
    fig.savefig('desra.png')
    plt.close(fig)

    # fig,ax = plt.subplots()
    # ax.set_xlabel('step')
    # ax.set_ylabel('gamma', labelpad=4.0)
    # ax.plot(np.array(model.gamma_list)*100,lw=1) #label
    # plt.show()
    # plt.close(fig)

    # fig,ax = plt.subplots()
    # ax.set_xlabel('cycle')
    # ax.set_ylabel('p', labelpad=4.0)
    # ax.plot(cycle,p,lw=1) #label
    # ax.set_ylim(0,p.max()*1.1)
    # plt.show()
    # plt.close(fig)

    # fig,ax = plt.subplots()
    # ax.set_xlabel('gamma')
    # ax.set_ylabel('tau', labelpad=4.0)
    # ax.plot(gamma*100,tau/1000,lw=1) #label
    # plt.show()
    # plt.close(fig)