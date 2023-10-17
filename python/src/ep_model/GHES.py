import numpy as np
import scipy.linalg
import math
import ep_model.GHES_hd_FEM as hd_FEM
import ep_model.GHES_hd3d_FEM as hd3d_FEM

class GHES:
    def __init__(self,G0,nu,phi=30.0,hmax=0.20):
        # Elastic parameters
        self.G0,self.nu = G0,nu
        # Parameters
        self.phi,self.hmax = phi,hmax
        self.K0 = 2*G0*(1+nu)/(3*(1-2*nu))

        self.p_ref = 40.e3 
        self.G = G0

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
        self.qmodel = hd3d_FEM.Multi_spring_3d(self.p_ref,G0,self.g05,hmax,model2d=hd_FEM.GHE_S)

        # identity
        self.Z3 = np.zeros([3,3])
        self.I3 = np.eye(3) 

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
    def isotropic_compression(self,e0,compression_stress):
        self.stress = self.I3 * compression_stress
        self.strain = np.zeros((3,3))

    # -------------------------------------------------------------------------------------- #
    def plastic_deformation(self,dstrain_given,dstress_given,deformation,sp):
        d = self.matrix_to_vector_bool(deformation)
        Ep = self.ep_modulus(sp)

        dstrain = self.matrix_to_vector(dstrain_given)
        dstrain[d] = 0.0
        dstress_constrain = np.dot(Ep,dstrain)
        dstress = self.matrix_to_vector(dstress_given) - dstress_constrain

        dstress_mask = dstress[d]
        Ep_mask = Ep[d][:,d]
        dstrain_mask = np.linalg.solve(Ep_mask,dstress_mask)

        dstrain[d] = dstrain_mask
        dstress = np.dot(Ep,dstrain)

        # update inner state #
        strain_vec = self.matrix_to_vector(sp.strain) + dstrain
        _ = self.qmodel.shear(strain_vec,sp.p)

        dstrain_mat = self.vector_to_matrix(dstrain)
        dstress_mat = self.vector_to_matrix(dstress)
        sp2 = self.StateParameters(sp.strain,sp.stress,dstrain_mat,dstress_mat)

        return dstrain_mat,dstress_mat,sp2

    # -------------------------------------------------------------------------------------- #
    def ep_modulus(self,sp,de=1.e-6):
        Ep = np.zeros([6,6])

        for i in range(6):
            dstrain = np.zeros(6)
            dstrain[i] = de
            strain_vec = self.matrix_to_vector(sp.strain) + dstrain
            dev_stress_vec = self.qmodel.shear_(strain_vec,sp.p)
            pstress_vec = self.pmodel(dstrain,sp.p)
            stress_vec = pstress_vec + dev_stress_vec

            dstress = stress_vec - self.matrix_to_vector(sp.stress)
            Ep[i,:] = dstress/de

        return Ep

    # -------------------------------------------------------------------------------------- #
    def pmodel(self,dstrain_vec,p):
        dev = dstrain_vec[0] + dstrain_vec[1] + dstrain_vec[2]
        dp = self.K0 * dev

        dstress_vec = np.zeros(6)
        dstress_vec[0] = dp + p
        dstress_vec[1] = dp + p
        dstress_vec[2] = dp + p

        return dstress_vec

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
