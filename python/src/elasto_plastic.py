from errno import EINPROGRESS
import numpy as np
import matplotlib.pyplot as plt

from ep_model import Li
from ep_model.DL1d import constitution as DL1d

import sys
import warnings
import traceback
warnings.filterwarnings('error')

def instantiateEP(dof,style,param):
    if style == 'ep_Li':
        return EP(dof,style,param)
    elif style == 'ep_DL1d':
        return EP_DL1d(dof,style,param)
    elif style == 'ep_DL1d_light':
        return EP_DL1d_Light(dof,style,param)
    elif style == 'ep_DL1d_ghe':
        return EP_DL1d_GHE(dof,style,param)
    else:
        print('Element style does not exist.')
        return


class EP:
    def __init__(self,dof,style,param):
        self.dof = dof
        self.style = style

        self.set_model(param)
        self.clear_strain()

    # -------------------------------------------------------------------------------------- #
    def set_model(self,param):
        if self.style == "ep_Li":
            nu,G0,M,self.e0,eg,d1,cohesion = param
            self.model = Li.Li2002(G0=G0,nu=nu,M=M,eg=eg,d1=d1,cohesion=cohesion)

        deformation_vec = np.array([False,False,False,False,False,False],dtype=bool)
        self.deformation = self.vector_to_matrix(deformation_vec)

    # -------------------------------------------------------------------------------------- #
    # convert FEM stress vector to stress matrix
    # Note: positive in compression!
    def FEMstress_to_matrix(self,FEMstress):
        if self.dof == 1:
            stress_vec = np.array([0.0,0.0,0.0,-FEMstress[0],-FEMstress[1],0.0])
        elif self.dof == 2:
            stress_vec = np.array([-FEMstress[0],-FEMstress[0],-FEMstress[1],0.0,0.0,-FEMstress[2]])
        elif self.dof == 3:
            stress_vec = np.array([-FEMstress[0],-FEMstress[0],-FEMstress[1], \
                                         -FEMstress[3], -FEMstress[4], -FEMstress[2]])

        return self.vector_to_matrix(stress_vec)

    def matrix_to_FEMstress(self,stress):
        stress_vec = self.matrix_to_vector(stress)

        if self.dof == 1:
            FEMstress = np.array([-stress_vec[3],-stress_vec[4]])
        elif self.dof == 2:
            FEMstress = np.array([-stress_vec[0],-stress_vec[2],-stress_vec[5]])
        elif self.dof == 3:
            FEMstress = np.array([-stress_vec[0],-stress_vec[2],-stress_vec[5],-stress_vec[3],-stress_vec[4]])

        return FEMstress

    # -------------------------------------------------------------------------------------- #
    # convert FEM strain vector to strain matrix
    # Note: positive in compression!
    def FEMstrain_to_matrix(self,FEMstrain):
        if self.dof == 1:
            strain_vec = np.array([0.0,0.0,0.0,-0.5*FEMstrain[0],-0.5*FEMstrain[1],0.0])
        elif self.dof == 2:
            strain_vec = np.array([-FEMstrain[0],0.0,-FEMstrain[1],0.0,0.0,-0.5*FEMstrain[2]])
        elif self.dof == 3:
            strain_vec = np.array([-FEMstrain[0],0.0,-FEMstrain[1], \
                                         -0.5*FEMstrain[3], -0.5*FEMstrain[4], -0.5*FEMstrain[2]])

        return self.vector_to_matrix(strain_vec)

    def matrix_to_FEMstrain(self,strain):
        strain_vec = self.matrix_to_vector(strain)

        if self.dof == 1:
            FEMstrain = np.array([-2*strain_vec[3],-2*strain_vec[4]])
        elif self.dof == 2:
            FEMstrain = np.array([-strain_vec[0],-strain_vec[2],-2*strain_vec[5]])
        elif self.dof == 3:
            FEMstrain = np.array([-strain_vec[0],-strain_vec[2],-2*strain_vec[5],-2*strain_vec[3],-2*strain_vec[4]])

        return FEMstrain

    # -------------------------------------------------------------------------------------- #
    # convert material modulus (3,3,3,3) to D matrix (ndof,ndof)
    def modulus_to_Dmatrix(self,E):
        if self.dof == 1:
            D = np.zeros([2,2])
            D[0,0],D[0,1] = E[0,1,0,1],E[0,1,1,2]
            D[1,0],D[1,1] = E[1,2,0,1],E[1,2,1,2]

        elif self.dof == 2:
            D = np.zeros([3,3])
            D[0,0],D[0,1],D[0,2] = E[0,0,0,0],E[0,0,2,2], 0.5*(E[0,0,2,0]+E[0,0,0,2])
            D[1,0],D[1,1],D[1,2] = E[2,2,0,0],E[2,2,2,2], 0.5*(E[2,2,2,0]+E[2,2,0,2])
            D[2,0] = 0.5*(E[2,0,0,0]+E[0,2,0,0])
            D[2,1] = 0.5*(E[2,0,2,2]+E[0,2,2,2])
            D[2,2] = 0.25*(E[2,0,2,0]+E[0,2,2,0]+E[2,0,0,2]+E[0,2,0,2])

        elif self.dof == 3:
            D = np.zeros([5,5])
            D[0,0],D[0,1],D[0,2],D[0,3],D[0,4] = E[0,0,0,0],E[0,0,2,2],E[0,0,2,0],E[0,0,0,1],E[0,0,1,2]
            D[1,0],D[1,1],D[1,2],D[1,3],D[1,4] = E[2,2,0,0],E[2,2,2,2],E[2,2,2,0],E[2,2,0,1],E[2,2,1,2]
            D[2,0],D[2,1],D[2,2],D[2,3],D[2,4] = E[2,0,0,0],E[2,0,2,2],E[2,0,2,0],E[2,0,0,1],E[2,0,1,2]
            D[3,0],D[3,1],D[3,2],D[3,3],D[3,4] = E[0,1,0,0],E[0,1,2,2],E[0,1,2,0],E[0,1,0,1],E[0,1,1,2]
            D[4,0],D[4,1],D[4,2],D[4,3],D[4,4] = E[1,2,0,0],E[1,2,2,2],E[1,2,2,0],E[1,2,0,1],E[1,2,1,2]

        return D

    # -------------------------------------------------------------------------------------- #
    def vector_to_matrix(self,vec):
        mat = np.array([[vec[0],vec[3],vec[5]],
                        [vec[3],vec[1],vec[4]],
                        [vec[5],vec[4],vec[2]]])
        return mat

    def matrix_to_vector(self,mat):
        vec = np.array([mat[0,0],mat[1,1],mat[2,2],0.5*(mat[0,1]+mat[1,0]),0.5*(mat[1,2]+mat[2,1]),0.5*(mat[2,0]+mat[0,2])])
        return vec

    def clear_strain(self):
        self.strain = np.zeros((3,3))

    # -------------------------------------------------------------------------------------- #
    def elastic_modulus(self):
        p,_ = self.model.set_stress_variable(self.stress)
        G0,K0 = self.model.elastic_modulus(self.e,p)
        return G0, K0 - G0*2/3

    def elastic_modulus_ep(self,e=None,p=None):
        if e is None:
            e = self.e0
        G0,K0 = self.model.elastic_modulus(e,p)
        return G0, K0 - G0*2/3


    # -------------------------------------------------------------------------------------- #
    def initial_state(self,init_stress,last=None):
        init_stress_mat = self.FEMstress_to_matrix(init_stress)

        # isotropic_compression
        compression_stress = init_stress_mat[0,0]
        self.model.isotropic_compression(self.e0,compression_stress)
        self.model.e0 = np.copy(self.e0)
        self.model.e = np.copy(self.e0)

        # set initial parameters
        self.stress = self.model.stress
        self.strain = self.model.strain
        self.e = self.model.e

        p,_ = self.model.set_stress_variable(self.stress)
        self.model.beta,self.model.H2 = p,p

        nstep = 50
        dstrain_vec = np.zeros(6,dtype=np.float64)
        dstress_vec = np.zeros(6,dtype=np.float64)

        dstress_vec[2] = (init_stress_mat[2,2]-init_stress_mat[0,0])/nstep
        deformation_vec = np.array([True,True,True,True,True,True],dtype=bool)

        dstrain_input = self.vector_to_matrix(dstrain_vec)
        dstress_input = self.vector_to_matrix(dstress_vec)
        deformation = self.vector_to_matrix(deformation_vec)

        sp0 = self.model.StateParameters(self.strain,self.stress,dstrain_input,dstress_input,self.model.stress_shift)

        # load to triaxial stress state
        for i in range(nstep):
            sp = self.model.StateParameters(self.strain,self.stress,sp0.dstrain,dstress_input,self.model.stress_shift)

            p,R = self.model.set_stress_variable(self.stress)
            dstrain,dstress,sp0 = \
                self.model.plastic_deformation(dstrain_input,dstress_input,deformation,sp)

            self.stress += dstress
            self.strain += dstrain

            # print(self.stress)

            ev,gamma = self.model.set_strain_variable(self.strain)
            self.e = self.e0 - ev*(1+self.e0)
            # print(p,R,gamma)

        # FEMstress = self.matrix_to_FEMstress(self.stress)
        # FEMstrain = self.matrix_to_FEMstrain(self.strain)
        # print("strain: ",FEMstrain)
        # print("stress: ",FEMstress)
        # print("strain: ",self.strain.diagonal())
        # print("stress: ",self.stress.diagonal())

    # -------------------------------------------------------------------------------------- #
    def set_Dp_matrix(self,FEMdstrain):
        dstrain = self.FEMstrain_to_matrix(FEMdstrain)

        # print(" ")

        dstress_vec = np.zeros(6,dtype=np.float64)
        dstress_input = self.vector_to_matrix(dstress_vec)

        sp0 = self.model.StateParameters(self.strain,self.stress,dstrain,dstress_input,self.model.stress_shift)
        ef1,ef2 = self.model.check_unload(sp0)

        # print(ef1,self.model.alpha)
        # print(ef2,sp0.p,self.model.beta)


        sp = self.model.StateParameters(self.strain,self.stress,dstrain,dstress_input,self.model.stress_shift,ef1=ef1,ef2=ef2)
        Ep = self.model.plastic_stiffness(sp)

        dstrain,dstress = \
            self.model.solve_strain_with_consttain(dstrain,dstress_input,Ep,self.deformation)

        sp2 = self.model.StateParameters(self.strain,self.stress,dstrain,dstress,self.model.stress_shift,ef1=ef1,ef2=ef2)
        self.model.update_parameters(sp2)

        ev,gamma = self.model.set_strain_variable(self.strain)
        self.e = self.e0 - ev*(1+self.e0)

        self.strain += dstrain
        self.stress += dstress

        # print(self.strain)
        # print(self.stress)

        Dp = self.modulus_to_Dmatrix(Ep)
        stress_yy = -self.stress[1,1]

        return Dp, self.matrix_to_FEMstress(self.stress), stress_yy

    # -------------------------------------------------------------------------------------- #
    def strain_to_stress(self,FEMdstrain):
        dstrain = self.FEMstrain_to_matrix(FEMdstrain)

        dstress_vec = np.zeros(6,dtype=np.float64)
        dstress_input = self.vector_to_matrix(dstress_vec)

        sp = self.model.StateParameters(self.strain,self.stress,dstrain,dstress_input)
        p,R = self.model.set_stress_variable(self.stress)

        dstrain,dstress,_ = \
            self.model.plastic_deformation(dstrain,dstress_input,self.deformation,sp)

        ev,gamma = self.model.set_strain_variable(self.strain)
        self.e = self.e0 - ev*(1+self.e0)

        self.strain += dstrain
        self.stress += dstress

        # p_stress,_ = np.linalg.eig(self.stress)
        # print(np.min(p_stress),np.max(p_stress))

        return self.matrix_to_FEMstress(self.stress)


class EP_DL1d(EP):
    '''
    Args:
        `style (string)`: 'ep_DL1d'
        `param (tuple)`: Tuple of `(nu,rho,info)`
        `info (dict)`: Dictionary that contains material parameters. Keys are 'H', 'P0', 'N', 'G0', 'sand', 'silt', 'clay', 'wL' and 'wP'. The first four are in SI units and the last five are in percent(%).
    '''
    def __init__(self,dof,style,param):
        super().__init__(dof,style,param)

    def set_model(self,param):
        self.rho,self.nu,self.info = param
        self.G0 = self.info['G0']
        self.G = self.G0
        self.K = 2*self.G*(1+self.nu)/3/(1-2*self.nu)
        self.p0 = self.info['P0']
        self.model = DL1d.DL1d(self.info)

    def elastic_modulus(self,G=None):
        if G is None:
            G = self.G
        K = self.K
        rlambda = K-G*2/3
        # rlambda = 2*G*self.nu/(1+self.nu)
        return G,rlambda

    def elastic_modulus_ep(self, e=None, p=None):
        return self.elastic_modulus()

    def initial_state(self, init_stress, last=False):
        init_stress_mat = self.FEMstress_to_matrix(init_stress)
        init_stress_vec = self.matrix_to_vector(init_stress_mat)

        p = init_stress_vec[:3].mean() * 1e-3  # N -> kN
        print(f'\t\tp/p0:{p/self.p0}')
        self.G = self.G0 * np.sqrt(p/self.p0)
        self.K = 2*self.G*(1+self.nu)/3/(1-2*self.nu)
        # self.G = self.G0

        E_inv = 1/(2*self.G*(1+self.nu))
        nu_by_E = -self.nu*E_inv
        mu_inv = 1/(2*self.G)
        S_mat = np.array([
            [E_inv,nu_by_E,nu_by_E,0,0,0],
            [nu_by_E,E_inv,nu_by_E,0,0,0],
            [nu_by_E,nu_by_E,E_inv,0,0,0],
            [0,0,0,mu_inv,0,0],
            [0,0,0,0,mu_inv,0],
            [0,0,0,0,0,mu_inv]
        ])

        init_strain_vec = S_mat@init_stress_vec
        self.strain = self.vector_to_matrix(init_strain_vec)
        self.stress = init_stress_mat
        self.gamma_list,self.tau_list = [0],[0]
        self.elastic_flag = True
        self.elastic_limit = 1e-7

        if last:
            h = p*1e3 / (9.8*self.rho)  # elementから直接与える
            h *= 1.5
            info_update = {'G0':self.G,'P0':p,'H':h}
            # info_update = {}
            self.model.initial_state(info_update)

    def set_Dp_matrix(self, FEMdstrain):
        def get_Dp(rmu,rlambda):
            if self.dof == 1:
                print('1-dof not set yet')
            elif self.dof == 2:
                Dp = np.array([
                    [rlambda+2*rmu,rlambda,0],
                    [rlambda,rlambda+2*rmu,0],
                    [0,0,rmu]
                ])
            elif self.dof == 3:
                Dp = np.array([
                    [rlambda+2*rmu,rlambda,rlambda,0,0,0],
                    [rlambda,rlambda+2*rmu,rlambda,0,0,0],
                    [rlambda,rlambda,rlambda+2*rmu,0,0,0],
                    [0,0,0,rmu,0,0],
                    [0,0,0,0,rmu,0],
                    [0,0,0,0,0,rmu]
                ])
            Dp_half = np.array([
                [rlambda+2*rmu,rlambda,rlambda,0,0,0],
                [rlambda,rlambda+2*rmu,rlambda,0,0,0],
                [rlambda,rlambda,rlambda+2*rmu,0,0,0],
                [0,0,0,2*rmu,0,0],
                [0,0,0,0,2*rmu,0],
                [0,0,0,0,0,2*rmu]
            ])
            return Dp,Dp_half

        dstrain = self.FEMstrain_to_matrix(FEMdstrain)
        dgamma = 2*dstrain[0,2]
        if np.abs(dgamma)>1e-10:
            dtau,gamma = self.model.shear_d(dgamma)
            if self.elastic_flag:
                if np.abs(gamma)<self.elastic_limit:
                    self.G = self.G0
                else:
                    self.G = np.min([np.abs(dtau/dgamma),self.G0])
                    self.elastic_flag = False
            else:
                self.G = np.min([np.abs(dtau/dgamma),self.G0])
        else:
            dtau = 0.0
        self.gamma_list.append(self.gamma_list[-1]+dgamma)
        self.tau_list.append(self.tau_list[-1]+dtau)
        rmu,rlambda = self.elastic_modulus(self.G)
        Dp,Dp_half = get_Dp(rmu,rlambda)
        dstrain_vec = self.matrix_to_vector(dstrain)
        dstress_vec = Dp_half@dstrain_vec

        self.strain += dstrain
        self.stress += self.vector_to_matrix(dstress_vec)

        stress_yy = -self.stress[1,1]

        return Dp,self.matrix_to_FEMstress(self.stress),stress_yy

    def plot(self,fname='result/constitution_ep.png'):
        gamma = 100*np.array(self.gamma_list)
        tau = 1/1000*np.array(self.tau_list)
        fig,ax = plt.subplots(figsize=[5,5])
        ax.set_xlabel('gamma')
        ax.set_ylabel('tau', labelpad=4.0)
        ax.plot(gamma,tau) #label
        fig.savefig(fname)
        plt.close(fig)


class EP_DL1d_Light(EP_DL1d):
    def set_model(self,param):
        super().set_model(param)
        self.model = DL1d.DL1d(self.info,maxlen=1000)

    def plot(self,fname='result/constitution_ep_light.png'):
        super().plot(fname)


class EP_DL1d_GHE(EP_DL1d):
    def set_model(self,param):
        super().set_model(param)
        self.model = DL1d.GHE_mix(self.info)

    def plot(self,fname='result/constitution_ep_ghe.png'):
        super().plot(fname)





# --------------------------------#
if __name__ == "__main__":
    dof = 2
    param = [0.33,202,0.97,0.6975,0.957,0.0] # nu, G0, M, e0, eg, d1
    ep = EP(dof,"ep_Li",param)


    K0 = 0.5
    rho = 1700.0
    z = 5.0

    sx = -rho*9.8*z*K0
    sz = -rho*9.8*z
    init_stress = np.array([sx,sz,0.0])

    ep.initial_state(init_stress)

    dg = 0.00001
    strain = np.zeros(3)

    gamma_list = []
    tau_list = []
    for i in range(1000):
        dstrain = np.array([0.0,0.0,dg])
        strain += dstrain
        stress = ep.strain_to_stress(dstrain)
        print(strain[2],stress)

        gamma_list += [strain[2]]
        tau_list += [stress[2]]


    plt.figure()
    plt.plot(gamma_list,tau_list)
    plt.show()
