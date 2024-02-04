import numpy as np
import scipy.linalg
import math, copy, sys, os
import ep_model.GHES_hd_FEM as hd_FEM
import ep_model.GHES_hd3d_FEM as hd3d_FEM


class GHES:
    def __init__(self,G0,nu,phi=30.0,hmax=0.20,p_ref=40.e3):
        # Elastic parameters
        self.G0,self.nu = G0,nu
        # Parameters
        self.phi,self.hmax = phi,hmax
        self.K0 = 2*G0*(1+nu)/(3*(1-2*nu))

        self.p_ref = p_ref
        self.G = G0
        self.K = self.K0

        # stress parameters
        self.pr = 101.e3
        self.pmin = 1.e-4
        self.stress_shift = 0.0

        # strength parameter
        phi = np.deg2rad(self.phi)
        self.Mc = 6*np.sin(phi)/(3-np.sin(phi))
        tauf = self.Mc*self.p_ref
        self.gr = tauf/self.G0
        self.g05 = self.gr/2.5

        self.e0 = 0.80
        self.e = 0.80

        # model 
        self.qmodel = hd3d_FEM.Multi_spring_3d(self.p_ref,G0,self.g05,hmax,model2d=hd_FEM.GHE_S,ntheta=3)
        # self.qmodel = hd3d_FEM.Multi_spring_3d(self.p_ref,G0,self.g05,hmax,model2d=hd_FEM.GHE,ntheta=16)

        # identity
        self.Z3 = np.zeros([3,3])
        self.I3 = np.eye(3) 
        self.epsilon = sys.float_info.epsilon

        # # check file
        # id = 0
        # file_name = "tmp/ghes_"+'{0:02}'.format(id)
        # if os.path.isfile(file_name):
        #     os.remove(file_name)
        # self.file = open(file_name,"w") 

    # def __del__(self):
    #     self.file.close()

    # -------------------------------------------------------------------------------------- #
    class StateParameters:
        def __init__(self,strain,stress,dstrain,dstress,stress_shift=0.0,pore_pressure=0.0,ef1=False,ef2=False):
            self.strain = np.copy(strain)
            self.stress = np.copy(stress)
            self.dstress = np.copy(dstress)
            self.dstrain = np.copy(dstrain)
            self.pore_pressure = pore_pressure
            self.pmin = 1.e-4

            self.set_stress_variable(stress_shift)
            self.set_stress_increment()

            self.elastic_flag1 = ef1
            self.elastic_flag2 = ef2

        def set_stress_variable(self,stress_shift):
            self.p = np.trace(self.stress)/3
            self.sij = self.stress - self.p*np.eye(3)
            self.rij = self.sij / max(self.p+stress_shift,self.pmin)
            self.R = np.sqrt(1.5*np.square(self.rij).sum())

        def set_stress_increment(self):
            stress = self.stress + self.dstress
            p = np.trace(stress)/3
            self.dp = p - self.p


    # -------------------------------------------------------------------------------------- #
    def elastic_modulus(self,e,p):
        # G = self.G0
        G = self.G0/self.p_ref * p 
        K = G*2*(1+self.nu)/(3*(1-2*self.nu))  
        return G,K

    # -------------------------------------------------------------------------------------- #
    def set_strain_variable(self,strain):
        ev = np.trace(strain)
        dev_strain = strain - ev/3.0 * self.I3
        gamma = math.sqrt(2.0/3.0)*np.linalg.norm(dev_strain)
        return ev,gamma

    # -------------------------------------------------------------------------------------- #
    def set_stress_variable(self,stress):
        p = np.trace(stress)/3
        r_stress = (stress - p*self.I3) / max(p,self.pmin)
        R = math.sqrt(1.5)*np.linalg.norm(r_stress)
        return p,R

    # -------------------------------------------------------------------------------------- #
    def check_unload(self,sp):
        return False,False

    # -------------------------------------------------------------------------------------- #
    def plastic_stiffness(self,sp):
        dstrain_vec = self.matrix_to_vector(sp.dstrain)
        strain_update_vec = self.matrix_to_vector(sp.strain) + dstrain_vec
        stress_update_vec = self.strain_to_stress_(strain_update_vec)
        dstress_residual = stress_update_vec - self.matrix_to_vector(sp.stress)

        # np.set_printoptions(precision=5)
        # print("dstress_residual",dstress_residual)

        Ep = self.ep_modulus(sp,dstrain_vec)
        dstress_trial = np.dot(Ep,dstrain_vec)

        # print("dstress_trial",dstress_trial)

        for i in range(6):
            # if np.abs(dstress_trial[i]) > 1.e-8:
            if np.abs(dstress_trial[i]) > 1.e-6:
                Ep[i,:] = Ep[i,:] * dstress_residual[i]/dstress_trial[i]

        # dstress_trial = np.dot(Ep,dstrain_vec)
        # print("dstress_trial",dstress_trial)
        # print("")

        return Ep

    # -------------------------------------------------------------------------------------- #
    def solve_strain_with_constrain(self,strain_given,stress_given,E,deformation):
        # deformation: True => deform (stress given), False => constrain (strain given)
        d = deformation.flatten()
        A = np.reshape(E,(9,9))

        strain = np.copy(strain_given.flatten())
        strain[d] = 0.0                        # [0.0,0.0,given,...]
        stress_constrain = np.dot(A,strain)
        stress = np.copy(stress_given.flatten()) - stress_constrain

        stress_mask = stress[d]
        A_mask = A[d][:,d]

        Ainv = np.linalg.pinv(A_mask)
        strain_mask = Ainv @ stress_mask

        strain[d] = strain_mask
        stress = np.dot(A,strain)
        return np.reshape(strain,(3,3)), np.reshape(stress,(3,3))

    # -------------------------------------------------------------------------------------- #
    def update_parameters(self,sp):
        strain_vec = self.matrix_to_vector(sp.strain + sp.dstrain)
        _ = self.qmodel.shear(strain_vec,sp.p)

    # -------------------------------------------------------------------------------------- #
    def isotropic_compression(self,e0,compression_stress):
        self.stress = self.I3 * compression_stress
        self.strain = np.zeros((3,3))
        self.pc = compression_stress

    # -------------------------------------------------------------------------------------- #
    def plastic_deformation(self,dstrain_given,dstress_given,deformation,sp0):
        d = self.matrix_to_vector_bool(deformation)
        
        strain_target = sp0.strain + dstrain_given
        stress_target = sp0.stress + dstress_given
        dstrain_residual = np.copy(dstrain_given)
        dstress_residual = np.copy(dstress_given)
        # print("stress target",stress_target)

        sp = self.StateParameters(sp0.strain,sp0.stress,dstrain_residual,dstress_residual)

        for it in range(10):
            # dstress = Ep * dstrain #
            Ep = self.ep_modulus(sp)

            dstrain = self.matrix_to_vector(dstrain_residual)
            dstrain[d] = 0.0
            dstress_constrain = np.dot(Ep,dstrain)
            dstress = self.matrix_to_vector(dstress_residual) - dstress_constrain

            dstress_mask = dstress[d]
            Ep_mask = Ep[d][:,d]
            
            dstrain_mask = np.linalg.solve(Ep_mask,dstress_mask)

            dstrain[d] = dstrain_mask
            dstress = np.dot(Ep,dstrain)

            # check stress+ = Ep(strain+) 
            dstrain_mat = self.vector_to_matrix(dstrain)
            dstress_mat = self.vector_to_matrix(dstress)
            strain_update = sp.strain + dstrain_mat
            strain_update_vec = self.matrix_to_vector(strain_update)

            stress_update_vec = self.strain_to_stress_(strain_update_vec)

            stress_update = self.vector_to_matrix(stress_update_vec)

            dstrain_residual = strain_target - strain_update
            dstress_residual = stress_target - stress_update

            # print("stress update",stress_update_vec)
            sp = self.StateParameters(strain_update,stress_update,dstrain_residual,dstress_residual)

            norm = np.linalg.norm(dstress_residual) 
            if norm < 1.e-6:
                break

        # update inner state #
        stress_update_vec = self.strain_to_stress(strain_update_vec)

        dstrain_mat = strain_update - sp0.strain
        dstress_mat = stress_update - sp0.stress
        sp = self.StateParameters(sp0.strain,sp0.stress,dstrain_mat,dstress_mat)

        return dstrain_mat,dstress_mat,sp

    # -------------------------------------------------------------------------------------- #
    def update_parameters(self,sp):
        strain_update_vec = self.matrix_to_vector(sp.strain + sp.dstrain)
        self.strain_to_stress(strain_update_vec)
        
    # -------------------------------------------------------------------------------------- #
    def modulus_to_Dmatrix(self,E):
        D = np.zeros([3,3])

        D[0,0],D[0,1],D[0,2] = E[0,0],E[0,2],E[0,5] 
        D[1,0],D[1,1],D[1,2] = E[2,0],E[2,2],E[2,5] 
        D[2,0],D[2,1],D[2,2] = E[5,0],E[5,2],E[5,5] 

        return D

    # -------------------------------------------------------------------------------------- #
    def ep_tensor(self,sp,de=1.e-6):
        Ep = np.zeros([3,3,3,3])

        for i in range(3):
            for j in range(3):
                dstrain = np.zeros([3,3])
                dstrain[i,j] = de
                dstrain_vec = self.matrix_to_vector(dstrain)
                strain_vec = self.matrix_to_vector(sp.strain + dstrain)
                dev_stress_vec = self.qmodel_shear_(strain_vec,sp.p)
                pstress_vec = self.pmodel(dstrain_vec,sp.p)
                stress_vec = pstress_vec + dev_stress_vec

                stress = self.vector_to_matrix(stress_vec)
                dstress = stress - sp.stress
                Ep[:,:,i,j] = dstress/de

        return Ep

    # -------------------------------------------------------------------------------------- #
    def ep_modulus(self,sp,de=None):
        Ep = self.e_modulus(sp.p)

        if de is None:
            de = np.ones(6) * 1.e-8

        for i in range(6):
            # if abs(de[i]) > self.epsilon:
            if abs(de[i]) > 1.e-9:
                dstrain = np.zeros(6)
                dstrain[i] = de[i]
                strain_vec = self.matrix_to_vector(sp.strain) + dstrain
                stress_vec = self.strain_to_stress_(strain_vec)
                dstress = stress_vec - self.matrix_to_vector(sp.stress)
                Ep[:,i] = dstress[:]/de[i]

        # for i in range(0,3):
        #     if abs(de[i]) > self.epsilon:
        #         dstrain = np.zeros(6)
        #         dstrain[i] = de[i]
        #         strain_vec = self.matrix_to_vector(sp.strain) + dstrain
        #         stress_vec = self.strain_to_stress_(strain_vec)
        #         dstress = stress_vec - self.matrix_to_vector(sp.stress)
        #         Ep[0:3,i] = dstress[0:3]/de[i]

        # for i in range(3,6):
        #     if abs(de[i]) > self.epsilon:
        #         dstrain = np.zeros(6)
        #         dstrain[i] = de[i]
        #         strain_vec = self.matrix_to_vector(sp.strain) + dstrain
        #         stress_vec = self.strain_to_stress_(strain_vec)
        #         dstress = stress_vec - self.matrix_to_vector(sp.stress)
        #         Ep[3:6,i] = dstress[3:6]/de[i]

        return Ep

    # -------------------------------------------------------------------------------------- #
    def e_modulus(self,p):
        E = np.zeros([6,6])
        G0,K0 = self.elastic_modulus(self.e,p)
        rmu, rlambda = G0, K0-G0*2/3

        E[0,0],E[0,1],E[0,2] = rlambda+2*rmu, rlambda, rlambda
        E[1,0],E[1,1],E[1,2] = rlambda, rlambda+2*rmu, rlambda
        E[2,0],E[2,1],E[2,2] = rlambda, rlambda, rlambda+2*rmu

        E[3,3] = rmu
        E[4,4] = rmu
        E[5,5] = rmu

        return E

    # -------------------------------------------------------------------------------------- #
    def strain_to_stress(self,strain_vec):
        pstress_vec = self.pmodel(strain_vec)
        p = pstress_vec[0]
        dev_stress_vec = self.qmodel_shear(strain_vec,p)
        stress_vec = pstress_vec + dev_stress_vec

        # output_line = "{} {}\n".format(strain_vec[0],stress_vec[0])
        # self.file.write(output_line)

        return stress_vec

    def strain_to_stress_(self,strain_vec):
        pstress_vec = self.pmodel(strain_vec)
        p = pstress_vec[0]
        dev_stress_vec = self.qmodel_shear_(strain_vec,p)
        stress_vec = pstress_vec + dev_stress_vec
        return stress_vec

    def strain_to_stress_check(self,strain_vec):
        pstress_vec = self.pmodel(strain_vec)
        p = pstress_vec[0]
        dev_stress_vec = self.qmodel_shear_(strain_vec,p)

        print("")
        np.set_printoptions(precision=12)
        print(strain_vec)
        np.set_printoptions(precision=5)
        print(p,dev_stress_vec)
        print("")

        stress_vec = pstress_vec + dev_stress_vec
        return stress_vec

    # -------------------------------------------------------------------------------------- #
    def qmodel_shear(self,strain_vec,p):
        dev_stress_vec = self.qmodel.shear(strain_vec,p)
        return dev_stress_vec/self.p_ref * p

    def qmodel_shear_(self,strain_vec,p):
        dev_stress_vec = self.qmodel.shear_(strain_vec,p)
        return dev_stress_vec/self.p_ref * p

    # -------------------------------------------------------------------------------------- #
    def pmodel(self,strain_vec):
        ev = strain_vec[0] + strain_vec[1] + strain_vec[2]
        p = self.pc + self.K0 * ev

        stress_vec = np.zeros(6)
        stress_vec[0] = p
        stress_vec[1] = p
        stress_vec[2] = p

        return stress_vec

    # -------------------------------------------------------------------------------------- #
    def vector_to_matrix(self,vec):
        mat = np.array([[vec[0],vec[3],vec[5]],
                        [vec[3],vec[1],vec[4]],
                        [vec[5],vec[4],vec[2]]])
        return mat

    def matrix_to_vector(self,mat):
        vec = np.array([mat[0,0],mat[1,1],mat[2,2],0.5*(mat[0,1]+mat[1,0]),0.5*(mat[1,2]+mat[2,1]),0.5*(mat[2,0]+mat[0,2])])
        return vec

    def matrix_to_vector_bool(self,mat):
        vec = np.array([mat[0,0],mat[1,1],mat[2,2],mat[0,1],mat[1,2],mat[2,0]])
        return vec
