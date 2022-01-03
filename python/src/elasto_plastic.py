import numpy as np
import matplotlib.pyplot as plt

from ep_model import Li

import sys

class EP:
    def __init__(self,dof,style,param):
        self.dof = dof
        self.style = style

        self.set_model(param)
        self.clear_strain()

    # -------------------------------------------------------------------------------------- #
    def set_model(self,param):
        if self.style == "ep_Li":
            nu,G0,M,self.e0,eg,d1 = param
            self.model = Li.Li2002(G0=G0,nu=nu,M=M,eg=eg,d1=d1)

        deformation_vec = np.array([False,False,False,False,False,False],dtype=bool)
        self.deformation = self.vector_to_matrix(deformation_vec)

    # -------------------------------------------------------------------------------------- #
    # convert FEM stress vector to stress matrix
    # Note: positive in compression!
    def FEMstress_to_matrix(self,FEMstress):
        if self.dof == 1:
            stress_vec = np.array([0.0,0.0,0.0,FEMstress[0],FEMstress[1],0.0])
        elif self.dof == 2:
            stress_vec = np.array([-FEMstress[0],-FEMstress[0],-FEMstress[1],0.0,0.0,FEMstress[2]])
        elif self.dof == 3:
            stress_vec = np.array([-FEMstress[0],-FEMstress[0],-FEMstress[1], \
                                         FEMstress[3], FEMstress[4], FEMstress[2]])

        return self.vector_to_matrix(stress_vec)

    def matrix_to_FEMstress(self,stress):
        stress_vec = self.matrix_to_vector(stress)

        if self.dof == 1:
            FEMstress = np.array([stress_vec[3],stress_vec[4]])
        elif self.dof == 2:
            FEMstress = np.array([-stress_vec[0],-stress_vec[2],stress_vec[5]])
        elif self.dof == 3:
            FEMstress = np.array([-stress_vec[0],-stress_vec[2],stress_vec[5],stress_vec[3],stress_vec[4]])

        return FEMstress

    # -------------------------------------------------------------------------------------- #
    # convert FEM strain vector to strain matrix
    # Note: positive in compression!
    def FEMstrain_to_matrix(self,FEMstrain):
        if self.dof == 1:
            strain_vec = np.array([0.0,0.0,0.0,0.5*FEMstrain[0],0.5*FEMstrain[1],0.0])
        elif self.dof == 2:
            strain_vec = np.array([-FEMstrain[0],0.0,-FEMstrain[1],0.0,0.0,0.5*FEMstrain[2]])
        elif self.dof == 3:
            strain_vec = np.array([-FEMstrain[0],0.0,-FEMstrain[1], \
                                         0.5*FEMstrain[3], 0.5*FEMstrain[4], 0.5*FEMstrain[2]])

        return self.vector_to_matrix(strain_vec)

    def matrix_to_FEMstrain(self,strain):
        strain_vec = self.matrix_to_vector(strain)

        if self.dof == 1:
            FEMstrain = np.array([2*strain_vec[3],2*strain_vec[4]])
        elif self.dof == 2:
            FEMstrain = np.array([-strain_vec[0],-strain_vec[2],2*strain_vec[5]])
        elif self.dof == 3:
            FEMstrain = np.array([-strain_vec[0],-strain_vec[2],2*strain_vec[5],2*strain_vec[3],2*strain_vec[4]])

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
            D[0,0],D[0,1],D[0,2] = E[0,0,0,0],E[0,0,2,2],E[0,0,2,0]
            D[1,0],D[1,1],D[1,2] = E[2,2,0,0],E[2,2,2,2],E[2,2,2,0]
            D[2,0],D[2,1],D[2,2] = E[2,0,0,0],E[2,0,2,2],E[2,0,2,0]

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
        vec = np.array([mat[0,0],mat[1,1],mat[2,2],mat[0,1],mat[1,2],mat[2,0]])
        return vec

    def clear_strain(self):
        self.strain = np.zeros((3,3))

    # -------------------------------------------------------------------------------------- #
    def elastic_modulus(self):
        p,_ = self.model.set_stress_variable(self.stress)
        G0,K0 = self.model.elastic_modulus(self.e,p)
        return G0, K0 - G0*2/3

    def elastic_modulus_ep(self,e,p):
        G0,K0 = self.model.elastic_modulus(e,p)
        return G0, K0 - G0*2/3


    # -------------------------------------------------------------------------------------- #
    def initial_state(self,init_stress):
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

        sp0 = self.model.StateParameters(self.strain,self.stress,dstrain_input,dstress_input)

        # load to triaxial stress state
        for i in range(nstep):
            sp = self.model.StateParameters(self.strain,self.stress,sp0.dstrain,dstress_input)

            p,R = self.model.set_stress_variable(self.stress)
            dstrain,dstress,sp0 = \
                self.model.plastic_deformation(dstrain_input,dstress_input,deformation,sp)

            self.stress += dstress
            self.strain += dstrain

            ev,gamma = self.model.set_strain_variable(self.strain)
            self.e = self.e0 - ev*(1+self.e0)
            # print(p,R,gamma)

        # FEMstress = self.matrix_to_FEMstress(self.stress)
        # FEMstrain = self.matrix_to_FEMstrain(self.strain)
        # print("strain: ",FEMstrain)
        # print("stress: ",FEMstress)

    # -------------------------------------------------------------------------------------- #
    def set_Dp_matrix(self,FEMdstrain):
        dstrain = self.FEMstrain_to_matrix(FEMdstrain)

        dstress_vec = np.zeros(6,dtype=np.float64)
        dstress_input = self.vector_to_matrix(dstress_vec)

        sp0 = self.model.StateParameters(self.strain,self.stress,dstrain,dstress_input)
        ef1,ef2 = self.model.check_unload(sp0)

        sp = self.model.StateParameters(self.strain,self.stress,dstrain,dstress_input,ef1=ef1,ef2=ef2)
        Ep = self.model.plastic_stiffness(sp)

        dstrain,dstress = \
            self.model.solve_strain_with_consttain(dstrain,dstress_input,Ep,self.deformation)

        ev,gamma = self.model.set_strain_variable(self.strain)
        self.e = self.e0 - ev*(1+self.e0)

        self.strain += dstrain
        self.stress += dstress

        Dp = self.modulus_to_Dmatrix(Ep)
        return Dp, self.matrix_to_FEMstress(self.stress)

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
