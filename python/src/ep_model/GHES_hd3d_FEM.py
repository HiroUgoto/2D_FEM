import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
try:
    import GHES_hd_FEM as hd_FEM
except:
    from . import GHES_hd_FEM as hd_FEM

np.set_printoptions(precision=5,suppress=True)

class Multi_spring_2d:
    def __init__(self,p0,G0,gr05=6e-5,hmax=0.25,model2d=hd_FEM.GHE_S,ntheta=12,idxy="",ghe_args=hd_FEM.S_PARAM):
        self.ntheta = ntheta
        self.dtheta = np.pi/ntheta
        self.theta = np.array([i*self.dtheta for i in range(ntheta)])
        self.model = [model2d(G0,gr05,hmax=hmax,id=i,idxy=idxy,**ghe_args) for i in range(ntheta)]

    def shear(self,gamma_axis,gamma_theta,p):
        gamma = [-gamma_axis*np.sin(theta) + gamma_theta*np.cos(theta) for theta in self.theta]
        tau = np.array([self.model[i].shear1(gamma[i]) for i in range(self.ntheta)])
        d = self.dtheta / (np.pi/2)
        tau_axis = tau * np.sin(-self.theta) * d
        tau_theta = tau * np.cos(self.theta) * d
        return tau_axis.sum(),tau_theta.sum()

    def shear_(self,gamma_axis,gamma_theta,p):
        gamma = [-gamma_axis*np.sin(theta) + gamma_theta*np.cos(theta) for theta in self.theta]
        tau = np.array([self.model[i].shear_(gamma[i]) for i in range(self.ntheta)])
        d = self.dtheta / (np.pi/2)
        tau_axis = tau * np.sin(-self.theta) * d
        tau_theta = tau * np.cos(self.theta) * d
        return tau_axis.sum(),tau_theta.sum()

class Multi_spring_3d:
    def __init__(self,p0,G0,gr05=6e-5,hmax=0.25,model2d=hd_FEM.GHE_S,ntheta=12,ghe_args=hd_FEM.S_PARAM):
        self.tauf = G0*2.5*gr05
        args = p0,G0,gr05,hmax,model2d,ntheta
        self.xy = Multi_spring_2d(*args,idxy="xy",ghe_args=ghe_args)
        self.yz = Multi_spring_2d(*args,idxy="yz",ghe_args=ghe_args)
        self.zx = Multi_spring_2d(*args,idxy="zx",ghe_args=ghe_args)

    def strain_to_shear(self,v):
        xy = v[0]-v[1], 2*v[3] ##式(3.8)
        yz = v[1]-v[2], 2*v[4]
        zx = v[2]-v[0], 2*v[5]
        # v0_list = []
        # v0_list += [v[0]]
        return xy,yz,zx

    def shear_to_stress(self,xy,yz,zx,p,effective=True):
        stress = np.zeros(6)
        stress[0] = -(zx[0]-xy[0])*2/3  ## 式(3.9)
        stress[1] = -(xy[0]-yz[0])*2/3
        stress[2] = -(yz[0]-zx[0])*2/3
        if not effective:
            stress[:3] += p
        stress[3] = xy[1]
        stress[4] = yz[1]
        stress[5] = zx[1]
        return stress

    def shear_to_q(self,xy,yz,zx):
        phi1 = 4*(xy[0]**2 + yz[0]**2 + zx[0]**2)
        phi2 = xy[1]**2 + yz[1]**2 + zx[1]**2
        q = np.sqrt(phi1/2 + 3*phi2)
        return q

    def shear(self,strain,p,ret_vec=True,ret_q=False,effective=True):
        exy,eyz,ezx = self.strain_to_shear(strain)
        sxy = self.xy.shear(*exy,p)
        syz = self.yz.shear(*eyz,p)
        szx = self.zx.shear(*ezx,p)
        vec,q = None,None
        if ret_vec:
            vec = self.shear_to_stress(sxy,syz,szx,p,effective)
        if ret_q:
            q = self.shear_to_q(sxy,syz,szx)
            # if abs(q)>self.tauf*1.5:
            #     print('!!!',abs(q)*1e-3,self.tauf*1e-3)
        if ret_vec and ret_q:
            return vec,q
        elif ret_vec:
            return vec
        else:
            return q
    
    def shear_(self,strain,p,ret_vec=True,ret_q=False,effective=True):
        exy,eyz,ezx = self.strain_to_shear(strain)
        sxy = self.xy.shear_(*exy,p)
        syz = self.yz.shear_(*eyz,p)
        szx = self.zx.shear_(*ezx,p)
        vec,q = None,None
        if ret_vec:
            vec = self.shear_to_stress(sxy,syz,szx,p,effective)
        if ret_q:
            q = self.shear_to_q(sxy,syz,szx)
            # if abs(q)>self.tauf*1.5:
            #     print('!!!',abs(q)*1e-3,self.tauf*1e-3)
        if ret_vec and ret_q:
            return vec,q
        elif ret_vec:
            return vec
        else:
            return q


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import dilatancy

    p0 = 100e3
    ncycle = 20
    div = 100
    amp = 0.15/2

    ncycle = 20
    # div = 50
    # amp = 0.1
    # amp = 0.02

    cycle = np.linspace(0,ncycle,ncycle*div)
    ea = amp * np.sin(cycle*2*np.pi)
    ea[:int(len(ea))] *= np.linspace(0,1,int(len(ea)))
    gamma = np.sqrt(3)*ea
    gamma = 1.5*ea
    strain_vec1 = np.zeros([len(gamma),6])
    strain_vec2 = np.zeros([len(gamma),6])
    strain_vec1[:,5] = gamma/2
    strain_vec2[:,0] = -ea/2
    strain_vec2[:,1] = -ea/2
    strain_vec2[:,2] = ea
    # strain_vec2[:,0] = gamma/2
    # strain_vec2[:,0] = -gamma/2
    # strain_vec1 = strain_vec2.copy()

    # G0 = 80e6
    G0 = 50e6
    gr = 6e-5*2.5
    gr = 1e-5
    print(G0*(2.5*gr)*1e-3)
    # multi1 = Multi_spring_3d(p0,G0/np.sqrt(p0/1000),model2d=hd_FEM.special_GHE_S,ntheta=2)
    # multi2 = Multi_spring_3d(p0,G0/np.sqrt(p0/1000),model2d=hd_FEM.special_GHE_S,ntheta=2)
    multi1 = Multi_spring_3d(p0,G0,gr,model2d=hd_FEM.special_GHE_S,ntheta=2)
    multi2 = Multi_spring_3d(p0,G0,gr,model2d=hd_FEM.special_GHE_S,ntheta=2)
    single1 = hd_FEM.special_GHE_S(G0,gr,0.25,**hd_FEM.S_PARAM)
    single2 = hd_FEM.special_GHE_S(G0,gr,0.25,**hd_FEM.S_PARAM)
    # single2 = hd_FEM.GHE_S(G0,gr,0.25,**hd_FEM.S_PARAM)
    tau1 = np.zeros_like(gamma)
    tau2 = np.zeros_like(gamma)
    tau3 = np.zeros_like(gamma)
    tau4 = np.zeros_like(gamma)
    q = np.zeros_like(gamma)
    p_model1 = dilatancy.Martin_Finn_Seed3(p0,k2=1e-4,M=1.5)
    p_model2 = dilatancy.Martin_Finn_Seed3(p0,k2=1e-4,M=1.5)
    p_list1 = np.zeros_like(gamma)
    p_list2 = np.zeros_like(gamma)
    stress_vec = np.zeros(6)
    for i in range(len(gamma)):
        stress_vec1,q[i] = multi1.shear(strain_vec1[i],p0,True,True)
        tau1[i] = stress_vec1[5]
        # tau1[i] = stress_vec1[2]-stress_vec1[0]

        stress_vec2 = multi2.shear(strain_vec2[i],p0,True)
        # tau2[i] = (stress_vec2[2]-stress_vec2[0])/np.sqrt(3)
        tau2[i] = (stress_vec2[2]-stress_vec2[0])/2
        # tau2[i] = stress_vec2[2]-stress_vec2[0]

        # tau3[i] = single1.shear(ea[i])
        tau4[i] = single2.shear(gamma[i])

        p_list1[i] = p_model1.shear(strain_vec1[i],stress_vec1)
        p_list2[i] = p_model2.shear(strain_vec2[i],stress_vec2)

    plt.rcParams.update({'lines.linewidth':1})

    fig,ax = plt.subplots()
    ax.plot(ea*100,tau1/1000,label='Shear',ls='-') #label
    ax.plot(ea*100,tau2/1000,label='Triaxial',ls='-.') #label
    # ax.plot(ea*100,tau3/1000,label='Single ea',ls='-.') #label
    ax.plot(ea*100,tau4/1000,label='Single gamma',ls='-.') #label
    ax.legend()
    plt.show()
    plt.close(fig)

    # fig,ax = plt.subplots()
    # ax.plot(p_list1/1000,tau1/1000,label='Shear') #label
    # ax.plot(p_list2/1000,tau2/1000,label='Triaxial',ls='-.') #label
    # ax.set_xlim(0,p_list1.max()/1000*1.1)
    # ax.legend()
    # plt.show()
    # plt.close(fig)

    # fig,ax = plt.subplots()
    # ax.plot(cycle,p_list1/1000,label='Shear') #label
    # ax.plot(cycle,p_list2/1000,label='Triaxial',ls='-.') #label
    # ax.set_ylim(0,p_list1.max()/1000*1.1)
    # ax.legend()
    # plt.show()
    # plt.close(fig)

    # M = 1.25
    # tauf = p0*M/np.sqrt(3)
    # # gr,hmax = 1e-3,0.2
    # gr,hmax = 1e-4,0.2
    # g05 = gr/2.5
    # G0 = tauf/g05
    # print('tau_f =',tauf/1000)https://www.kyoto-u.ac.jp/ja/news/2023-04-13-0 
    # models = []
    # # models.append(hd_FEM.mod_GHE(G0=G0,gr=g05,hmax=hmax,**hd_FEM.S_PARAM))
    # models.append(hd_FEM.special_GHE_S(G0=G0,gr=g05,hmax=hmax,**hd_FEM.S_PARAM))
    # model2d = hd_FEM.Multi(models)

    # model = GHE_3d(model2d)
    # (gamma_i,tau_i),(gamma_s,tau_s),d = model.tortion(gamma)

    # plt.rcParams.update({'lines.linewidth':1})

    # # fig,ax = plt.subplots()
    # # ax.plot(gamma*100,label='input') #label
    # # ax.plot(gamma_s*100,label='gamma_s') #label
    # # ax.plot(gamma_i*100,label='gamma_i') #label
    # # ax.plot(d/np.pi) #label
    # # ax.legend()
    # # plt.show()
    # # plt.close(fig)

    # # fig,ax = plt.subplots()
    # # ax.set_ylabel('tau_s', labelpad=4.0)
    # # ax.plot(tau_s) #label
    # # plt.show()
    # # plt.close(fig)

    # # fig,ax = plt.subplots()
    # # ax.set_xlabel('gamma_i [%]')
    # # ax.set_ylabel('tau_i [kN]', labelpad=4.0)
    # # ax.plot(gamma_i*100,tau_i/1000) #label
    # # plt.show()
    # # plt.close(fig)

    # fig,ax = plt.subplots()
    # ax.set_xlabel('gamma_s [%]')
    # ax.set_ylabel('tau_s [kN]', labelpad=4.0)
    # ax.plot(gamma_s*100,tau_s/1000) #label
    # plt.show()
    # plt.close(fig)