import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
try:
    import GHES_hd_FEM as hd_FEM
except:
    from . import GHES_hd_FEM as hd_FEM

np.set_printoptions(precision=5,suppress=True)

class GHE_3d:
    def __init__(self,model2d=hd_FEM.myGHE_S(50e3,5e-4,0.2)):
        self.model = model2d
        self.strain_vec = np.zeros(6)
        self.sign = 1
        self.eig_vec = np.eye(3)
        self.gamma_i = 0.0
        self.gamma_s = 0.0
        self.tau_i,self.tau_s = 0.0,0.0
        self.dgamma_i = 1.0
        self.reverse = False

    def vector_to_matrix(self,vec):
        if vec.ndim >= 2:
            mat = np.zeros(list(vec.shape[:2])+[3,3])
            # mat = np.zeros((len(vec),3,3))
            mat[...,0,0] = vec[...,0]
            mat[...,1,1] = vec[...,1]
            mat[...,2,2] = vec[...,2]
            mat[...,0,1] = vec[...,3]
            mat[...,1,0] = vec[...,3]
            mat[...,1,2] = vec[...,4]
            mat[...,2,1] = vec[...,4]
            mat[...,0,2] = vec[...,5]
            mat[...,2,0] = vec[...,5]
        else:
            mat = np.array([[vec[0],vec[3],vec[5]],
                            [vec[3],vec[1],vec[4]],
                            [vec[5],vec[4],vec[2]]])
        return mat

    def dq(self,r1,r2):
        dq = (r1*r2.inv()).as_quat()
        d = 2*np.arccos(np.abs(dq[-1]))
        return d

    def d_rot(self,strain_vec):
        # strain_mat.shape: (3,3)

        def det_sign(A):
            eijk = np.zeros([3,3,3])
            eijk[0,1,2] = eijk[1,2,0] = eijk[2,0,1] = 1
            eijk[0,2,1] = eijk[2,1,0] = eijk[1,0,2] = -1
            det = np.einsum('...ijk,...i,...j,...k',eijk,A[...,0,:],A[...,1,:],A[...,2,:])
            return np.sign(det)

        # def continuous(vec,prev_vec):
        def continuous(val,vec,prev_vec):
            # vec.shape: (3,3)
            vecs = np.stack([vec.copy() for i in range(4)],axis=0)
            vecs[1,1] *= -1
            vecs[1,2] *= -1
            vecs[2,0] *= -1
            vecs[2,2] *= -1
            vecs[3,0] *= -1
            vecs[3,1] *= -1
            rs = Rotation.from_matrix(vecs)
            prev_r = Rotation.from_matrix(prev_vec)
            d = []
            for r in rs:
                d.append(self.dq(r,prev_r))
            d = np.array(d)
            vec = vecs[np.argmin(d)]
            rmin = rs[np.argmin(d)]
            # print(eig_vec,(rmin*prev_r.inv()).as_quat(),d.min()>np.pi/3,'\n')
            # print(val,d,d.min()>np.pi/3,'\n')
            return vec,d.min()

        strain_mat = self.vector_to_matrix(strain_vec)
        eig_val,eig_vec = np.linalg.eig(strain_mat)
        idx = eig_val.argsort()[::-1]
        eig_val = eig_val[idx]
        eig_vec = eig_vec[:,idx]
        eig_vec = eig_vec.T

        eig_vec[...,2] *= det_sign(eig_vec)[...,np.newaxis]
        eig_vec,d = continuous(eig_val,eig_vec,self.eig_vec)
        # eig_vec,d = continuous(eig_vec,self.eig_vec)
        self.eig_vec = eig_vec
        # print(strain_mat,'\n',eig_vec,d>np.pi/3,'\n')
        return d

    def get_gamma(self,strain):
        e = strain.copy()
        e[:3] -= strain[:3].mean()
        gamma = 2*np.sqrt((e**2).sum())  # gamma for sheer
        return gamma


    def shear(self,strain_vec):
        d = self.d_rot(strain_vec-self.strain_vec)
        reverse = d>np.pi/3
        self.d = d

        gamma_i = self.get_gamma(strain_vec)
        dgamma_i = gamma_i-self.gamma_i
        # if (dgamma_i*self.dgamma_i<0) and (not reverse):
        # if reverse:
        s_reverse = dgamma_i*self.dgamma_i < 0
        if s_reverse and not reverse:
            self.sign *= -1
        dgamma_s = self.sign * dgamma_i
        self.gamma_s += dgamma_s
        self.gamma_i = gamma_i
        self.dgamma_i = dgamma_i

        self.strain_vec = strain_vec
        self.tau_s = self.model.shear(self.gamma_s)
        self.tau_i = np.abs(self.tau_s)
        self.reverse = reverse
        return self.tau_i

    def test(self,strain_vecs):
        gamma_i,gamma_s,tau_i,tau_s = [],[],[],[]
        d = []
        for strain_vec in strain_vecs:
            self.shear(strain_vec)
            gamma_i.append(self.gamma_i)
            gamma_s.append(self.gamma_s)
            tau_i.append(self.tau_i)
            tau_s.append(self.tau_s)
            d.append(self.d)
        return (np.array(gamma_i),np.array(tau_i)),(np.array(gamma_s),np.array(tau_s)),np.array(d)

    def tortion(self,gammas):
        seq = len(gammas)
        strain_vecs = np.zeros([seq,6])
        strain_vecs[:,5] = gammas/2
        return self.test(strain_vecs)
        # (gamma_i,tau_i),(gamma_s,tau_s) = self.test(strain_vecs)

    def triaxial(self,ax_strains):
        seq = len(ax_strains)
        strain_vecs = np.zeros([seq,6])
        strain_vecs[:,2] = ax_strains
        strain_vecs[:,0] = strain_vecs[:,1] = -ax_strains/2
        return self.test(strain_vecs)


class Multi_spring_2d:
    def __init__(self,p0,G0,gr05=6e-5,hmax=0.25,model2d=hd_FEM.GHE_S,ntheta=12,ghe_args=hd_FEM.S_PARAM):
        self.ntheta = ntheta
        self.dtheta = np.pi/ntheta
        self.theta = np.array([i*self.dtheta for i in range(ntheta)])
        self.model = [model2d(G0,gr05,hmax=hmax,id=i,**ghe_args) for i in range(ntheta)]

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
        self.xy = Multi_spring_2d(*args,ghe_args=ghe_args)
        self.yz = Multi_spring_2d(*args,ghe_args=ghe_args)
        self.zx = Multi_spring_2d(*args,ghe_args=ghe_args)

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