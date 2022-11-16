import numpy as np
import math
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class LD2004:
    # Defalt parameters are for Toyoura sand (LD2004)
    def __init__(self,G0=125,nu=0.25,M=1.25,c=0.75,eg0=0.86,rlambdac=0.019,xi=0.7, \
                 d1=0.41,m=3.5,h1=2.3,h2=2.23,h3=2.2,n=1.1, \
                 d2=1,h4=3.5,a=1,b1=0.005,b2=2,b3=0.01, \
                 kg=0.59,delta=0.2,kh=0.2,d3=1.1,d4=9.5,d5=0.12,h5=0.32,kr=0.6,\
                 cohesion=0.0):
        # Elastic parameters
        self.G0,self.nu = G0,nu
        # Critical state parameters
        self.M,self.c,self.rlambdac,self.xi = M,c,rlambdac,xi
        # parameters associated with dr-mechanisms
        self.d1,self.m,self.h1,self.h2,self.h3,self.n = d1,m,h1,h2,h3,n
        # parameters associated with dp-mechanisms
        self.d2,self.h4 = d2,h4
        # parameters associated with L
        self.a,self.b1,self.b2,self.b3 = a,b1,b2,b3

        # parameters associated with third mechanism
        self.delta = delta
        self.eg0,self.kg,self.kh = eg0,kg,kh
        self.d3,self.d4,self.d5,self.h5,self.kr = d3,d4,d5,h5,kr

        # minimum epsillon
        self.eps = 1.e-6

        # stress parameters
        self.pr = 101.e3
        self.pmin = 100.0

        # stress & strain
        self.stress = np.zeros((3,3))
        self.strain = np.zeros((3,3))

        # cohesion parameters
        self.cohesion = cohesion
        self.stress_shift = self.cohesion/self.M    # stress shift applied in State Parameter

        # BS parameters
        self.alpha = np.zeros((3,3))
        self.beta = 0.0
        self.H1 = 0.0
        self.H2 = 0.0

        # Accumulated index
        self.L1 = 0.0

        # Third paramneters
        self.eg = self.eg0
        self.h_factor = 1.0
        self.L3 = 0.0

        # Fabric tensor
        if self.delta == 1.0:
            self.Finv = np.diag([2,2,0])
        else:
            self.Finv = np.diag([1/(1+self.delta),1/(1+self.delta),1/(1-self.delta)]) * (3+self.delta)

        # identity
        self.Z3 = np.zeros([3,3])
        self.I3 = np.eye(3)
        self.Dijkl = np.einsum('ij,kl->ijkl',self.I3,self.I3)
        self.Dikjl = np.einsum('ij,kl->ikjl',self.I3,self.I3)
        self.Diljk = np.einsum('ij,kl->ilkj',self.I3,self.I3)

        # parameters
        self.sqrt2_3 = math.sqrt(2/3)
        self.fn = 2*(1+self.nu)/(3*(1-2*self.nu))
        self.rlambda_coeff = 2*self.nu/(1-2*self.nu)
        self.G2_coeff = self.fn*self.h4 / (self.h4 + np.sqrt(2/3)*self.fn*self.d2) / self.fn
        self.g0 = self.c*(1+self.c) / (1+self.c**2)
        self.dg0 = (-self.c**2*(1-self.c)*(1+self.c)**2) / (1+self.c**2)**3

    # -------------------------------------------------------------------------------------- #
    class StateParameters:
        def __init__(self,strain,stress,dstrain,dstress,stress_shift,ef1=False,ef2=False,ef3=False):
            self.strain = np.copy(strain)
            self.stress = np.copy(stress)
            self.dstress = np.copy(dstress)
            self.dstrain = np.copy(dstrain)
            self.pmin = 1.0

            self.set_stress_variable(stress_shift)
            self.set_stress_increment(stress_shift)
            self.set_stress_increment_tensor(stress_shift)

            self.elastic_flag1 = ef1
            self.elastic_flag2 = ef2
            self.elastic_flag3 = ef3

        def set_stress_variable(self,stress_shift):
            self.p = np.trace(self.stress)/3
            self.sij = self.stress - self.p*np.eye(3)
            self.rij = self.sij / max(self.p+stress_shift,self.pmin)
            self.R = np.sqrt(1.5*np.square(self.rij).sum())

        def set_stress_increment(self,stress_shift):
            stress = self.stress + self.dstress
            p = np.trace(stress)/3
            self.dp = p - self.p

        def set_stress_increment_tensor(self,stress_shift):
            d  = np.einsum('ik,jl->ijkl',np.eye(3),np.eye(3))
            sd = np.einsum('ij,kl->ijkl',self.stress,np.eye(3))
            self.M = (3*(self.p+stress_shift)*d - sd) / (3*(self.p+stress_shift)**2)
            self.drij = np.einsum('ijkl,kl->ij',self.M,self.dstress)

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
        self.strain = self.Z3.copy()

    # -------------------------------------------------------------------------------------- #
    def set_strain_variable(self,strain):
        ev = np.trace(strain)
        dev_strain = strain - ev/3.0 * self.I3
        gamma = math.sqrt(2.0/3.0)*np.linalg.norm(dev_strain)
        return ev,gamma

    def set_strain_variable_ev(self,strain):
        ev = np.trace(strain)
        return ev

    def set_strain_increment(self,dstrain_mat):
        strain_mat = self.strain_mat + dstrain_mat
        ev0,gamma0 = self.set_strain_variable(self.strain_mat)
        ev,gamma = self.set_strain_variable(strain_mat)
        return ev-ev0,gamma-gamma0

    # -------------------------------------------------------------------------------------- #
    def set_stress_variable(self,stress):
        p = np.trace(stress)/3
        r_stress = (stress - p*self.I3) / max(p,self.pmin)
        R = math.sqrt(1.5)*np.linalg.norm(r_stress)
        return p,R

    def set_stress_variable_p(self,stress):
        p = np.trace(stress)/3
        return p

    # -------------------------------------------------------------------------------------- #
    def Lode_angle(self,dev_stress):
        J2 = 0.5*np.square(dev_stress).sum()
        if J2 == 0.0:
            return 0.0
        J3 = -np.linalg.det(dev_stress)
        s3 = J3/2.0 * (3.0/J2)**1.5
        s3 = max(s3,-1.0)
        s3 = min(s3,1.0)
        theta3 = np.arcsin(s3)
        return theta3

    def g_theta(self,dev_stress):            # Eq.(7)
        theta3 = self.Lode_angle(dev_stress)
        if theta3 == 0.0:
            return self.g0
        st = np.sin(theta3)
        g1 = np.sqrt((1+self.c**2)**2 + 4*self.c*(1-self.c**2)*st) - (1+self.c**2)
        g2 = 2*(1-self.c)*st
        return g1/g2

    def dg_theta(self,theta3):                 # Eq.(45)
        st = np.sin(theta3)
        if st == 0.0:
            return self.dg0
        g11 = self.c*(1+self.c)
        g12 = st*np.sqrt((1+self.c**2)**2+4*self.c*(1-self.c**2)*st)
        g21 = np.sqrt((1+self.c**2)**2 + 4*self.c*(1-self.c**2)*st) - (1+self.c**2)
        g22 = 2*(1-self.c)*st * st
        dg = g11/g12 - g21/g22
        return dg

    def Lode_angle_J2(self,dev_stress):
        J2 = 0.5*np.square(dev_stress).sum()
        if J2 == 0.0:
            return 0.0,0.0
        J3 = -np.linalg.det(dev_stress)
        s3 = J3/2 * (3/J2)**1.5
        s3 = max(s3,-1.0)
        s3 = min(s3,1.0)
        theta3 = np.arcsin(s3)
        return theta3,J2

    def g_theta_J2(self,dev_stress):            # Eq.(7)
        theta3,J2 = self.Lode_angle_J2(dev_stress)
        if theta3 == 0.0:
            return self.g0,J2
        st = np.sin(theta3)
        g1 = np.sqrt((1+self.c**2)**2 + 4*self.c*(1-self.c**2)*st) - (1+self.c**2)
        g2 = 2*(1-self.c)*st
        return g1/g2,J2

    # -------------------------------------------------------------------------------------- #
    def state_parameter(self,e,p):
        psi = e - (self.eg-self.rlambdac*(max(p,self.pmin)/self.pr)**self.xi)  # Eq.(18)
        return psi

    # -------------------------------------------------------------------------------------- #
    def elastic_modulus(self,e,p):
        G = self.G0*(2.97-e)**2 / (1+e) * np.sqrt(max(p,self.pmin)*self.pr)  # Eq.(16)
        K = G*2*(1+self.nu)/(3*(1-2*self.nu))                           # Eq.(17)
        return G,K

    def elastic_modulus_G(self,e,p):
        G = self.G0*(2.97-e)**2 / (1+e) * np.sqrt(max(p,self.pmin)*self.pr)  # Eq.(16)
        return G

    def elastic_stiffness(self,G):
        rlambda = G*self.rlambda_coeff
        Ee = rlambda*self.Dijkl + G*(self.Dikjl+self.Diljk)
        return Ee

    def isotropic_compression_stiffness(self,e,p):
        G = self.elastic_modulus_G(e,p)
        G2 = G*self.G2_coeff  #  Elastic + Eq.(29)
        E2 = self.elastic_stiffness(G2)
        return E2

    # -------------------------------------------------------------------------------------- #
    def set_mapping_stress(self,sp,elastic_flag1=False,elastic_flag2=False):
        def mapping_r(t,rij,alpha):
            rij_bar = alpha + t*(rij-alpha)
            g_bar,J2 = self.g_theta_J2(rij_bar)
            R_bar = math.sqrt(3.0*J2)
            return rij_bar,R_bar,g_bar

        def mapping_r_Rg(t,rij,alpha):
            rij_bar = alpha + t*(rij-alpha)
            g_bar,J2 = self.g_theta_J2(rij_bar)
            R_bar = math.sqrt(3.0*J2)
            return R_bar,g_bar

        def F1_boundary_surface(t,*args):               # Eq.(6)
            rij,alpha = args
            R_bar,g_bar = mapping_r_Rg(t,rij,alpha)
            return R_bar-self.H1*g_bar

        def F1_boundary_surface_all(t,*args):               # Eq.(6)
            rij,alpha = args
            rij_bar,R_bar,g_bar = mapping_r(t,rij,alpha)
            return R_bar-self.H1*g_bar,rij_bar,R_bar,g_bar

        def find_rho1_ratio(rij,alpha):
            if np.abs(F1_boundary_surface(1.0,rij,alpha)) < 1.e-8:
                return 1.0
            t0,t1 = 0.99,1.e2
            F1_t0 = F1_boundary_surface(t0,rij,alpha)
            F1_t1 = F1_boundary_surface(t1,rij,alpha)
            for itr in range(10):
                if F1_t0*F1_t1 >= 0.0:
                    t1 = 10*t1
                    F1_t1 = F1_boundary_surface(t1,rij,alpha)
                else:
                    break
            try:
                t = scipy.optimize.brenth(F1_boundary_surface,t0,t1,args=(rij,alpha))
            except ValueError:
                print("ValueError",t0,t1,F1_t0,F1_t1)
                exit()

            return t

        H1,H2 = self.H1,self.H2

        if np.linalg.norm(sp.rij-self.alpha) < 1.e-6:  # Elastic behavior
            sp.elastic_flag1 = True
        else:
            F1,rij_bar,R_bar,g_bar = F1_boundary_surface_all(1.0,sp.rij,self.alpha)
            if F1 > 0.0:
                sp.rij_bar,sp.R_bar,sp.g_bar = rij_bar,R_bar,g_bar
                H1 = sp.R_bar/sp.g_bar
                sp.rho1_ratio = 1.0
            else:
                t = find_rho1_ratio(sp.rij,self.alpha)
                sp.rij_bar,sp.R_bar,sp.g_bar = mapping_r(t,sp.rij,self.alpha)
                sp.rho1_ratio = np.copy(t)      # rho1_ratio = rho1_bar / rho1

        if np.abs(sp.p-self.beta) < 1.e-6:  # Elastic behavior
            sp.elastic_flag2 = True
        else:
            if sp.p > self.H2:
                H2 = np.copy(sp.p)

            if sp.p > self.beta:
                sp.p_bar = np.copy(self.H2)
            else:
                sp.p_bar = np.copy(self.pmin)

            rho2 = np.abs(sp.p-self.beta)
            rho2_b = np.abs(sp.p_bar-self.beta)
            sp.rho2_ratio = max(rho2_b/rho2,1.0)

        return H1, H2

    # -------------------------------------------------------------------------------------- #
    def modified_stress_tensor_hat(self,sp):
        if sp.elastic_flag1:
            st = self.I3
        else:
            st = self.M*sp.g_bar/sp.R_bar * sp.rij_bar + self.I3

        T_hat = (st@self.Finv + self.Finv@st) /6.0
        p_hat = np.trace(T_hat)/3
        rij_hat = (T_hat - p_hat*self.I3) / p_hat
        R_hat = math.sqrt(1.5)*np.linalg.norm(rij_hat)
        g_hat = self.g_theta(rij_hat)

        A_bar = R_hat/(self.M*g_hat) - 1.0
        self.eg = self.eg0 + self.kg * A_bar * np.abs(A_bar)

        Ac = R_hat / self.M - 1.0
        Ae = R_hat / (self.M * self.c) - 1.0
        self.h_factor = ((self.kh*Ac - Ae) + (1-self.kh)*A_bar) / (Ac-Ae)

    def modified_stress_tensor_tilde(self,sp):
        if sp.R == 0.0:
            st = self.I3
        else:
            st = self.M*sp.g/sp.R * sp.rij + self.I3

        T_tilde = (st@self.Finv + self.Finv@st) /6.0
        p_tilde = np.trace(T_tilde)/3
        rij_tilde = (T_tilde - p_tilde*self.I3) / p_tilde
        R_tilde = math.sqrt(1.5)*np.linalg.norm(rij_tilde)
        g_tilde = self.g_theta(rij_tilde)

        A = R_tilde/(self.M*g_tilde) - 1.0
        Ac = R_tilde / self.M - 1.0
        Ae = R_tilde / (self.M * self.c) - 1.0
        return ((self.kr*Ac - Ae) + (1-self.kr)*A) / (Ac-Ae)

    # -------------------------------------------------------------------------------------- #
    def set_parameters(self,sp):
        def accumulated_load_index(L1):         # Eq.(22)
            fL = (1-self.b3)/np.sqrt((1-L1/self.b1)**2+(L1/self.b1)/self.b2**2) + self.b3
            return fL

        def scaling_factor(e,rho1_ratio):         # Eq.(21)
            fL = accumulated_load_index(self.L1)
            r1 = (1.0/rho1_ratio)**10
            h = (self.h1-self.h2*self.e)*(r1+self.h3*fL*(1-r1)) * self.h_factor
            return h

        def plastic_modulus1(G,R_bar,g_bar,rho1_ratio,h,psi):    # Eq.(19)
            Mg_R = self.M*g_bar*np.exp(-self.n*psi) / R_bar
            Kp1 = G*h*(Mg_R*rho1_ratio - 1)
            Kp1_b = G*h*(Mg_R - 1)
            return Kp1,Kp1_b

        def dilatancy1(R,g,rho1_ratio,psi):      # Eq.(23)
            R_Mg = R / (self.M*g)
            D1 = self.d1*(np.exp(self.m*psi)*np.sqrt(rho1_ratio) - R_Mg)
            return D1

        def plastic_modulus2(G,Mg_R,rho2_ratio,sign):   #  Eq.(25)
            Kp2 = G*self.h4*Mg_R * (rho2_ratio)**self.a*sign
            if rho2_ratio == 1.0 and sign > 0.0:
                Kp2_b = np.copy(Kp2)
            else:
                Kp2_b = 0.0
            return Kp2,Kp2_b

        def dilatancy2(Mg_R,sign):               #  Eq.(27)
            if Mg_R >= 1.0:
                D2 = self.d2*(Mg_R-1.0)*sign
            else:
                D2 = 0.0
            return D2

        def plastic_modulus3(G,Mg_R,factor):     # Eq.(25)
            Kp3 = G*self.h5/(self.pr*self.delta) * (Mg_R-1) * factor
            return Kp3

        def dilatancy3(psi,g,R):    # Eq.(24)
            D3 = self.d3 * np.exp(self.d4*psi) * (self.M*g-R) / (np.exp(self.d5*self.L3))
            return D3

        sp.Ge,sp.Ke = self.elastic_modulus(self.e,sp.p)
        sp.psi = self.state_parameter(self.e,sp.p)
        sp.g = self.g_theta(sp.sij)

        if sp.elastic_flag1:
            sp.Kp1_b = 0.0
            sp.D1 = 0.0
        else:
            h = scaling_factor(self.e,sp.rho1_ratio)
            sp.Kp1,sp.Kp1_b = plastic_modulus1(sp.Ge,sp.R_bar,sp.g_bar,sp.rho1_ratio,h,sp.psi)
            sp.D1 = dilatancy1(sp.R,sp.g,sp.rho1_ratio,sp.psi)

        if sp.elastic_flag2 or sp.R == 0.0:
            sp.Kp2_b = 0.0
            sp.D2 = 0.0
        else:
            sign = (sp.p-self.beta)/np.abs(sp.p-self.beta)
            Mg_R = self.M*sp.g/sp.R
            sp.Kp2,sp.Kp2_b = plastic_modulus2(sp.Ge,Mg_R,sp.rho2_ratio,sign)
            sp.D2 = dilatancy2(Mg_R,sign)

        if sp.elastic_flag3 or sp.R == 0.0:
            sp.Kp3_inv = 0.0
            sp.D3 = 0.0
            sp.alpha1 = 1.0
            sp.alpha2 = 0.0
        else:
            gm = 1.0 - sp.R/(self.M * sp.g)     # Eq.(20)
            if gm < 0.0:
                gm = 0.0
            sp.alpha1 = gm / math.sqrt(1-2*gm*(1-gm))       # Eq.(19)
            sp.alpha2 = (1-gm) / math.sqrt(1-2*gm*(1-gm))   # Eq.(19)

            Mg_R = self.M*sp.g/sp.R
            factor = self.modified_stress_tensor_tilde(sp)
            sp.Kp3_inv = 1.0 / plastic_modulus3(sp.Ge,Mg_R,factor)
            sp.D3 = dilatancy3(sp.psi,sp.g,sp.R)

    # -------------------------------------------------------------------------------------- #
    def set_parameter_nm(self,sp):
        def dF1_r(r_bar,R_bar,theta3_bar,g_bar,dg_bar):
            if np.abs(R_bar) < self.eps:
                return self.Z3
            st_bar = np.sin(theta3_bar)
            a = R_bar*g_bar + 3*R_bar*st_bar*dg_bar
            b = 9*dg_bar
            c = 1.5/(R_bar*g_bar)**2
            rr_bar = r_bar @ r_bar.T
            return (a*r_bar + b*rr_bar)*c

        if sp.elastic_flag1:
            sp.nij = np.eye(3)
        else:
            theta3_bar = self.Lode_angle(sp.sij)
            dg_bar = self.dg_theta(theta3_bar)
            dF1 = dF1_r(sp.rij_bar,sp.R_bar,theta3_bar,sp.g_bar,dg_bar)
            dF1_tr = np.trace(dF1)
            nij = dF1 - self.I3*dF1_tr/3.0
            sp.nij = nij / np.linalg.norm(nij)

        r_abs = math.sqrt(np.square(sp.rij).sum())
        if r_abs == 0.0:
            sp.mij = self.Z3.copy()
        else:
            sp.mij = sp.rij / r_abs         # mij = rij/|rij|

        if sp.elastic_flag3 or sp.R == 0.0:
            sp.lij_dash = np.zeros([3,3])
            sp.lij = np.zeros([3,3])
        else:
            ndr = np.einsum("ij,ij",sp.nij,sp.drij)
            drij_dash = sp.drij - ndr * sp.nij
            drij_dash_norm = np.linalg.norm(drij_dash)
            if drij_dash_norm == 0.0:
                sp.lij_dash = np.zeros([3,3])
            else:
                sp.lij_dash = drij_dash / drij_dash_norm
            sp.lij = sp.alpha1*sp.lij_dash + sp.alpha2*sp.nij

    # -------------------------------------------------------------------------------------- #
    def set_parameter_TZ(self,sp):
        if sp.elastic_flag1 or sp.elastic_flag2 or sp.R == 0.0:
            sp.B = 0.0
        else:
            nm = np.einsum("ij,ij",sp.nij,sp.mij)
            nr = np.einsum("ij,ij",sp.nij,sp.rij)
            Bu = 2*sp.Ge*nm - self.sqrt2_3*sp.Ke*sp.D2*nr
            Bd = self.sqrt2_3*sp.Ke*sp.D2 + sp.Kp2
            sp.B = Bu / Bd

        if sp.elastic_flag1:
            sp.Tij = self.Z3.copy()
        else:
            nr = np.einsum("ij,ij",sp.nij,sp.rij)
            Tu = 2*sp.Ge*sp.nij - sp.Ke*(nr+sp.B)*self.I3
            Td = 2*sp.Ge - self.sqrt2_3*sp.Ke*sp.D1*(nr+sp.B) + sp.Kp1
            sp.Tij = Tu / Td

        if sp.elastic_flag2:
            sp.Zij = self.Z3.copy()
        else:
            Zu = sp.Ke*self.I3 - self.sqrt2_3*sp.Ke*sp.D1*sp.Tij
            if sp.R == 0.0:
                Kp2_D2 = sp.Ge*self.h4/self.d2 * sp.rho2_ratio**self.a
                Zd = self.sqrt2_3*sp.Ke + Kp2_D2
            else:
                Zd = self.sqrt2_3*sp.Ke*sp.D2 + sp.Kp2
            sp.Zij = Zu / Zd

        if sp.elastic_flag3 or sp.R == 0.0:
            sp.B2,sp.B3 = 0,0
        else:
            nr = np.einsum("ij,ij",sp.nij,sp.rij)
            nl = np.einsum("ij,ij",sp.nij,sp.lij)

            if sp.elastic_flag1:
                sp.B2 = 0.0
            else:
                B2u = self.sqrt2_3*sp.Ke*sp.D3*(nr+sp.B) - 2*sp.Ge*nl
                B2d = 2*sp.Ge - self.sqrt2_3*sp.Ke*sp.D1*(nr+sp.B) + sp.Kp1
                sp.B2 = B2u / B2d

            if sp.elastic_flag2:
                sp.B3 = 0.0
            else:
                B3u = self.sqrt2_3*sp.Ke*(sp.D1*sp.B2 + sp.D3)
                B3d = self.sqrt2_3*sp.Ke*sp.D2 + sp.Kp2
                sp.B3 = B3u / B3d

    # -------------------------------------------------------------------------------------- #
    def set_parameter_LP(self,sp):
        def tensor_inv(A,B):
            C = np.zeros([3,3,3,3])
            A_reshape = np.reshape(A,(9,9))
            A_reshape_inv = np.linalg.pinv(A_reshape)
            for i in range(3):
                for j in range(3):
                    B_flat = B[:,:,i,j].flatten()
                    C_flat = A_reshape_inv @ B_flat
                    C[:,:,i,j] = np.reshape(C_flat,(3,3))
            return C

        Lm0 = (np.einsum('kp,lq->klpq',self.I3,self.I3) + np.einsum('kq,lp->klpq',self.I3,self.I3)) / 2.0
        if sp.elastic_flag1:
            Lm1 = np.zeros([3,3,3,3])
            Lm1d = np.zeros([3,3])
        else:
            nD = sp.nij + np.sqrt(2/27)*sp.D1*self.I3
            Lm1 = np.einsum("kl,pq->klpq",nD,sp.Tij)
            Lm1d = nD * sp.B2

        if sp.elastic_flag2:
            Lm2 = np.zeros([3,3,3,3])
            Lm2d = np.zeros([3,3])
        elif sp.R == 0:
            mD = np.sqrt(2/27)*self.I3
            Lm2 = np.einsum("kl,pq->klpq",mD,sp.Zij)
            Lm2d = mD * sp.B3
        else:
            mD = sp.mij + np.sqrt(2/27)*sp.D2*self.I3
            Lm2 = np.einsum("kl,pq->klpq",mD,sp.Zij)
            Lm2d = mD * sp.B3

        if sp.elastic_flag3:
            Lm3d = np.zeros([3,3])
            P = np.zeros([3,3])
        else:
            Lm3d = sp.lij + np.sqrt(2/27)*sp.D3*self.I3
            ldM = np.einsum("ij,ijkl->kl",sp.lij_dash,sp.M)
            P = sp.alpha1*ldM * sp.Kp3_inv

        Lm = Lm0 - Lm1 - Lm2
        Lmd = Lm1d - Lm2d + Lm3d
        Ee = self.elastic_stiffness(sp.Ge)
        Lm_dash = np.einsum('ijkl,klpq->ijpq',Ee,Lm)
        Lm2_dash = np.einsum('ijkl,kl->ij',Ee,Lmd)

        LP = Lm0 + np.einsum('ij,kl->ijkl',Lm2_dash,P)
        sp.L = tensor_inv(LP,Lm_dash)
        sp.PLij = np.einsum('ij,ijkl->kl',P,sp.L)

    # -------------------------------------------------------------------------------------- #
    def set_tensor_Ep(self,sp):
        sp.Ep = np.copy(sp.L)

    # -------------------------------------------------------------------------------------- #
    def check_unload(self,sp):
        self.set_mapping_stress(sp,False,False)
        self.modified_stress_tensor_hat(sp)
        self.set_parameters(sp)
        self.set_parameter_nm(sp)
        self.set_parameter_TZ(sp)
        self.set_parameter_LP(sp)

        dL3 = np.einsum("ij,ij",sp.PLij,sp.dstrain)
        if sp.elastic_flag1 or (dL3 < 0.0):
            elastic_flag3 = True
        else:
            elastic_flag3 = False

        dL1 = np.einsum("ij,ij",sp.Tij,sp.dstrain) + sp.B2*dL3
        if sp.elastic_flag1 or (dL1 < 0.0):
            elastic_flag1 = True
            self.alpha = np.copy(sp.rij)
        else:
            elastic_flag1 = False

        dL2 = np.einsum("ij,ij",sp.Zij,sp.dstrain) - sp.B3*dL3
        if sp.elastic_flag2 or (dL2 < -1.e-12):
            elastic_flag2 = True
            self.beta = np.copy(sp.p)
        else:
            elastic_flag2 = False

        return elastic_flag1,elastic_flag2,elastic_flag3

    # -------------------------------------------------------------------------------------- #
    def update_parameters(self,sp):
        sp.stress += sp.dstress
        H1,H2 = self.set_mapping_stress(sp)
        self.modified_stress_tensor_hat(sp)
        self.set_parameters(sp)
        self.set_parameter_nm(sp)
        self.set_parameter_TZ(sp)
        self.set_parameter_LP(sp)

        dL3 = np.einsum("ij,ij",sp.PLij,sp.dstrain)
        dL1 = np.einsum("ij,ij",sp.Tij,sp.dstrain) + sp.B2*dL3
        dL2 = np.einsum("ij,ij",sp.Zij,sp.dstrain) - sp.B3*dL3
        # print(" dL:",sp.elastic_flag1,dL1,sp.elastic_flag2,dL2)

        if not sp.elastic_flag1:
            self.L1 += dL1
            self.H1 = H1

        if not sp.elastic_flag2:
            self.H2 = H2

        if not sp.elastic_flag3:
            self.L3 += dL3 / sp.Kp3_inv

    # -------------------------------------------------------------------------------------- #
    def plastic_stiffness(self,sp):
        self.set_mapping_stress(sp)
        self.modified_stress_tensor_hat(sp)
        self.set_parameters(sp)
        self.set_parameter_nm(sp)
        self.set_parameter_TZ(sp)
        self.set_parameter_LP(sp)
        self.set_tensor_Ep(sp)
        return sp.Ep

    # -------------------------------------------------------------------------------------- #
    def solve_strain(self,stress_mat,E):
        b = stress_mat.flatten()
        A = np.reshape(E,(9,9))
        Ainv = np.linalg.pinv(A)
        x = Ainv @ b
        strain_mat = np.reshape(x,(3,3))
        return strain_mat

    def solve_strain_with_consttain(self,strain_given,stress_given,E,deformation):
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
    def elastic_deformation(self,dstrain,dstress,Ge,deformation):
        Ee = self.elastic_stiffness(Ge)
        dstrain_elastic,dstress_elastic = self.solve_strain_with_consttain(dstrain,dstress,Ee,deformation)
        return dstrain_elastic,dstress_elastic

    def plastic_deformation(self,dstrain_given,dstress_given,deformation,sp0):
        ef1,ef2,ef3 = self.check_unload(sp0)
        strain,stress = sp0.strain,sp0.stress

        sp = self.StateParameters(strain,stress,dstrain_given,dstress_given,self.stress_shift,ef1=ef1,ef2=ef2,ef3=ef3)
        Ep = self.plastic_stiffness(sp)
        dstrain_ep,dstress_ep = self.solve_strain_with_consttain(dstrain_given,dstress_given,Ep,deformation)

        sp2 = self.StateParameters(strain,stress,dstrain_ep,dstress_ep,self.stress_shift,ef1=ef1,ef2=ef2,ef3=ef3)
        self.update_parameters(sp2)

        dstrain = np.copy(dstrain_ep)
        dstress = np.copy(dstress_ep)

        return dstrain,dstress,sp2

    # -------------------------------------------------------------------------------------- #
    def isotropic_compression(self,e0,compression_stress,nstep=1000):
        dcp = compression_stress / nstep
        self.e = np.copy(e0)

        dstress_vec = np.array([dcp,dcp,dcp,0,0,0])
        dstress = self.vector_to_matrix(dstress_vec)

        self.stress = np.zeros((3,3))
        self.strain = np.zeros((3,3))
        for i in range(nstep):
            p = self.set_stress_variable_p(self.stress)
            E = self.isotropic_compression_stiffness(self.e,p)
            dstrain = self.solve_strain(dstress,E)

            self.stress += dstress
            self.strain += dstrain
            ev = self.set_strain_variable_ev(self.strain)
            self.e = e0 - ev*(1+e0)

        self.clear_strain()

    # -------------------------------------------------------------------------------------- #
    def triaxial_compression(self,e0,compression_stress,de=0.0001,emax=0.20,print_result=False,plot=False):
        self.isotropic_compression(e0,compression_stress)
        self.e0 = np.copy(e0)
        self.e = np.copy(e0)

        p,_ = self.set_stress_variable(self.stress)
        self.beta,self.H2 = p,p
        G0,K0 = self.elastic_modulus(self.e,p)

        nstep = int(emax/de)
        dstrain_vec = np.array([0.0,0.0,de,0.0,0.0,0.0])
        dstrain_input = self.vector_to_matrix(dstrain_vec)

        dstress_vec = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        dstress_input = self.vector_to_matrix(dstress_vec)

        deformation_vec = np.array([True,True,False,True,True,True],dtype=bool)
        deformation = self.vector_to_matrix(deformation_vec)

        sp0 = self.StateParameters(self.strain,self.stress,dstrain_input,dstress_input,self.stress_shift)

        gamma_list,R_list = [],[]
        ev_list = []
        max_p_stress = 0.0
        min_p_stress = compression_stress
        for i in range(nstep):
            sp = self.StateParameters(self.strain,self.stress,sp0.dstrain,dstress_input,self.stress_shift)

            p,R = self.set_stress_variable(self.stress)
            dstrain,dstress,sp0 = \
                self.plastic_deformation(dstrain_input,dstress_input,deformation,sp)

            self.stress += dstress
            self.strain += dstrain

            ev,gamma = self.set_strain_variable(self.strain)
            self.e = self.e0 - ev*(1+self.e0)

            p_stress,_ = np.linalg.eig(self.stress)

            # print(gamma,R,ev,p)
            # print(gamma,p_stress[0],p_stress[2])
            min_p_stress = min(min_p_stress,p_stress[0])
            max_p_stress = max(max_p_stress,p_stress[2])

            gamma_list += [gamma]
            R_list += [R]
            ev_list += [ev]

        if print_result:
            print("+++ triaxial_compression +++")
            print(" e0:",self.e0, "  e:",self.e)
            print(" G0[MPa]:",G0*1e-6)
            print(" sigma1:",max_p_stress, " sigma3:",min_p_stress)

        if plot:
            plt.figure()
            plt.plot(gamma_list,R_list)
            plt.show()

            plt.plot(gamma_list,ev_list)
            plt.show()

        return max_p_stress, min_p_stress, gamma_list, ev_list

    # -------------------------------------------------------------------------------------- #
    def cyclic_shear_test(self,e0,compression_stress,sr=0.2,cycle=10,print_result=False,plot=False):
        def cycle_load(nstep,amp,i):
            tau = amp*np.sin((i+1)/nstep*2*np.pi)
            tau0 = amp*np.sin((i+0)/nstep*2*np.pi)
            return tau - tau0

        self.isotropic_compression(e0,compression_stress)
        self.e0 = np.copy(e0)
        self.e = np.copy(e0)

        p0,_ = self.set_stress_variable(self.stress)
        self.beta,self.H2 = p0,p0

        nstep = 1000
        ncycle = cycle

        dstrain_vec = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        dstrain_input = self.vector_to_matrix(dstrain_vec)

        dstress_vec = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        dstress_input = self.vector_to_matrix(dstress_vec)

        deformation_vec = np.array([False,False,False,True,True,True],dtype=bool)
        deformation = self.vector_to_matrix(deformation_vec)

        sp0 = self.StateParameters(self.strain,self.stress,dstrain_input,dstress_input,self.stress_shift)

        gamma_list,tau_list = [],[]
        ev_list,p_list = [],[]
        step_list,ep_list = [],[]
        for ic in range(ncycle):
            print("N :",ic+1)
            for i in range(nstep):
                dtau = cycle_load(nstep,sr*p0,i)
                dstress_vec = np.array([0.0,0.0,0.0,0.0,0.0,dtau])
                dstress_input = self.vector_to_matrix(dstress_vec)

                sp = self.StateParameters(self.strain,self.stress,sp0.dstrain,dstress_input,self.stress_shift)

                p,R = self.set_stress_variable(self.stress)
                dstrain,dstress,sp0 = \
                    self.plastic_deformation(dstrain_input,dstress_input,deformation,sp)

                self.stress += dstress
                self.strain += dstrain

                ev,gamma = self.set_strain_variable(self.strain)
                self.e = self.e0 - ev*(1+self.e0)

                gamma_list += [self.strain[0,2]]
                tau_list += [self.stress[0,2]]
                p_list += [p]
                step_list += [ic*nstep+i]
                ep_list += [(p0-p)/p0]

        if plot:
            plt.figure()
            plt.plot(gamma_list,tau_list)
            plt.show()
            plt.plot(p_list,tau_list)
            plt.show()
            # plt.plot(step_list,ep_list)
            # plt.show()


    # -------------------------------------------------------------------------------------- #
    def cyclic_pure_shear_test(self,e0,compression_stress,gmax=0.01,cycle=10,print_result=False,plot=False):
        def cycle_load(nstep,amp,i):
            tau = amp*np.sin((i+1)/nstep*2*np.pi)
            tau0 = amp*np.sin((i+0)/nstep*2*np.pi)
            return tau - tau0

        self.isotropic_compression(e0,compression_stress*0.5)
        self.e0 = np.copy(e0)
        self.e = np.copy(e0)

        ###
        p,_ = self.set_stress_variable(self.stress)
        self.beta,self.H2 = p,p
        G0,K0 = self.elastic_modulus(self.e,p)

        nstep = 100
        ds = 0.5*compression_stress / nstep

        dstrain_vec = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        dstrain_input = self.vector_to_matrix(dstrain_vec)

        dstress_vec = np.array([0.0,0.0,ds,0.0,0.0,0.0])
        dstress_input = self.vector_to_matrix(dstress_vec)

        deformation_vec = np.array([True,True,True,True,True,True],dtype=bool)
        deformation = self.vector_to_matrix(deformation_vec)

        sp0 = self.StateParameters(self.strain,self.stress,dstrain_input,dstress_input,self.stress_shift)

        gamma_list,R_list = [],[]
        ev_list = []
        max_p_stress = 0.0
        min_p_stress = compression_stress
        for i in range(nstep):
            sp = self.StateParameters(self.strain,self.stress,sp0.dstrain,dstress_input,self.stress_shift)

            p,R = self.set_stress_variable(self.stress)
            dstrain,dstress,sp0 = \
                self.plastic_deformation(dstrain_input,dstress_input,deformation,sp)

            self.stress += dstress
            self.strain += dstrain

            ev,gamma = self.set_strain_variable(self.strain)
            self.e = self.e0 - ev*(1+self.e0)

            p_stress,_ = np.linalg.eig(self.stress)

            # print(gamma,R,ev,p)
            # print(gamma,p_stress[0],p_stress[2])
            min_p_stress = min(min_p_stress,p_stress[0])
            max_p_stress = max(max_p_stress,p_stress[2])

            gamma_list += [gamma]
            R_list += [R]
            ev_list += [ev]

        # print(self.strain)
        print(self.stress)
        self.clear_strain()
        ###

        nstep = 100
        ncycle = cycle

        dstrain_vec = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        dstrain_input = self.vector_to_matrix(dstrain_vec)

        dstress_vec = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        dstress_input = self.vector_to_matrix(dstress_vec)

        deformation_vec = np.array([False,False,True,False,False,False],dtype=bool)
        deformation = self.vector_to_matrix(deformation_vec)

        sp0 = self.StateParameters(self.strain,self.stress,dstrain_input,dstress_input,self.stress_shift)


        gamma_list,tau_list = [],[]
        ev_list,p_list = [],[]
        step_list,ep_list = [],[]
        for ic in range(ncycle):
            print("N :",ic+1)
            for i in range(nstep):
                dg = cycle_load(nstep,gmax,i)
                dstrain_vec = np.array([0.0,0.0,0.0,0.0,0.0,dg])
                dstrain_input = self.vector_to_matrix(dstrain_vec)

                sp = self.StateParameters(self.strain,self.stress,dstrain_input,sp0.dstress,self.stress_shift)

                p,R = self.set_stress_variable(self.stress)
                dstrain,dstress,sp0 = \
                    self.plastic_deformation(dstrain_input,dstress_input,deformation,sp)

                self.stress += dstress
                self.strain += dstrain

                ev,gamma = self.set_strain_variable(self.strain)
                self.e = self.e0 - ev*(1+self.e0)

                gamma_list += [self.strain[0,2]]
                tau_list += [self.stress[0,2]]
                p_list += [p]
                step_list += [ic*nstep+i]

                # print(self.strain)
            print(self.stress[0,0],self.stress[1,1],self.stress[2,2],self.stress[0,2])

        if plot:
            plt.figure()
            plt.plot(gamma_list,tau_list)
            plt.show()
            plt.plot(p_list,tau_list)
            plt.show()
            # plt.plot(step_list,ep_list)
            # plt.show()


# --------------------------------#
if __name__ == "__main__":

    emax = 0.957
    emin = 0.611
    Dr = 0.70
    e0 = emax-Dr*(emax-emin)
    print(e0)

    LD_model = LD2004(G0=420,nu=0.33,M=0.97,d1=0.41,delta=0.2,cohesion=4.e3)
    compression_stress = 40.e3
    LD_model.cyclic_pure_shear_test(e0,compression_stress,gmax=0.0001,cycle=10,print_result=True,plot=True)
    # LD_model.cyclic_pure_shear_test(e0,compression_stress,gmax=0.001,cycle=20,print_result=True,plot=True)
    exit()
