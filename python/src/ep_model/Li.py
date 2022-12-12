import numpy as np
import math
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class Li2002:
    # Defalt parameters are for Toyoura sand (Li2002)
    def __init__(self,G0=125,nu=0.25,M=1.25,c=0.75,eg=0.934,rlambdac=0.019,xi=0.7, \
                 d1=0.41,m=3.5,h1=2.1,h2=2.03,h3=2.2,n=1.1, \
                 d2=1.0,h4=3.5,a=1,b1=0.005,b2=2.0,b3=0.01,cohesion=0.0,e0=False):

        #parameter backup
        with open("./result/param.txt","w",newline="\n") as f:
            saveparam = "Li param\nG0={} nu={} \nM={} c={} eg={} rlambdac={} xi={} \nd1={} m={} h1={} h2={} h3={} n={} \nd2={} h4={} \na={} b1={} b2={} b3={} cohesion={} e0={}\n".format(G0,nu,M,c,eg,rlambdac,xi,d1,m,h1,h2,h3,n,d2,h4,a,b1,b2,b3,cohesion,e0)
            f.write(saveparam)
            f.close()

        # Elastic parameters
        self.G0,self.nu = G0,nu
        # Critical state parameters
        self.M,self.c,self.eg,self.rlambdac,self.xi = M,c,eg,rlambdac,xi
        # parameters associated with dr-mechanisms
        self.d1,self.m,self.h1,self.h2,self.h3,self.n = d1,m,h1,h2,h3,n
        # parameters associated with dp-mechanisms
        self.d2,self.h4 = d2,h4
        # Default parameters
        self.a,self.b1,self.b2,self.b3 = a,b1,b2,b3

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
        # self.G2_coeff = self.fn*self.h4 / (self.h4 + np.sqrt(2/3)*self.fn*self.d2) / self.fn
        self.G2_coeff = self.fn*self.h4 / (self.h4 + np.sqrt(2/3)*self.fn*self.d2)

        self.g0 = self.c*(1+self.c) / (1+self.c**2)
        self.dg0 = (-self.c**2*(1-self.c)*(1+self.c)**2) / (1+self.c**2)**3

        #CUtest
        self.epstress = 0.0
        self.fL = 1.0
        self.h = 0.0
        self.psi = 0.0
    # -------------------------------------------------------------------------------------- #
    class StateParameters:
        def __init__(self,strain,stress,dstrain,dstress,stress_shift,pore_pressure=0.0,ef1=False,ef2=False):
            self.strain = np.copy(strain)
            self.stress = np.copy(stress)
            self.dstress = np.copy(dstress)
            self.dstrain = np.copy(dstrain)
            self.pore_pressure = pore_pressure
            self.pmin = 1.0

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
        psi = e - (self.eg-self.rlambdac*(max(p+self.stress_shift,self.pmin)/self.pr)**self.xi)  # Eq.(18)
        self.psi = psi
        return psi

    # -------------------------------------------------------------------------------------- #
    def elastic_modulus(self,e,p):
        G = self.G0*(2.97-e)**2 / (1+e) * np.sqrt(max(p+self.stress_shift,self.pmin)*self.pr)  # Eq.(16)
        K = G*2*(1+self.nu)/(3*(1-2*self.nu))                           # Eq.(17)
        return G,K

    def elastic_modulus_G(self,e,p):
        G = self.G0*(2.97-e)**2 / (1+e) * np.sqrt(max(p+self.stress_shift,self.pmin)*self.pr)  # Eq.(16)
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
                # sp.p_bar = np.copy(self.pmin)
                sp.p_bar = np.copy(self.pmin-self.stress_shift)

            rho2 = np.abs(sp.p-self.beta)
            rho2_b = np.abs(sp.p_bar-self.beta)
            sp.rho2_ratio = max(rho2_b/rho2,1.0)

        return H1, H2

    # -------------------------------------------------------------------------------------- #
    def set_parameters(self,sp):
        def accumulated_load_index(L1):         # Eq.(22)
            fL = (1-self.b3)/np.sqrt((1-L1/self.b1)**2+(L1/self.b1)/self.b2**2) + self.b3
            self.fL = fL
            return fL

        def scaling_factor(e,rho1_ratio):         # Eq.(21)
            fL = accumulated_load_index(self.L1)
            r1 = (1.0/rho1_ratio)**10
            h = (self.h1-self.h2*self.e)*(r1+self.h3*fL*(1-r1))
            self.h = h
            return h

        def plastic_modulus1(G,R_bar,g_bar,rho1_ratio,h,psi):    # Eq.(19)
            Mg_R = self.M*g_bar*np.exp(-self.n*psi) / R_bar
            Kp1 = G*h*(Mg_R*rho1_ratio - 1)
            Kp1_b = G*h*(Mg_R - 1)
            return Kp1,Kp1_b

        def dilatancy1(R,g,rho1_ratio,psi):      # Eq.(23)
            R_Mg = R / (self.M*g)
            D1 = self.d1*(np.exp(self.m*psi)*np.sqrt(rho1_ratio) - R_Mg)
            # print("ep",np.exp(self.m*psi)*np.sqrt(rho1_ratio),"R_Mg",R_Mg)
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

        sp.Ge,sp.Ke = self.elastic_modulus(self.e,sp.p)
        sp.psi = self.state_parameter(self.e,sp.p)
        sp.g = self.g_theta(sp.sij)

        if sp.elastic_flag1:
            sp.Kp1_b = 0.0
            sp.D1 = 0.0
        else:
            h = scaling_factor(self.e,sp.rho1_ratio)
            sp.h = h
            sp.Kp1,sp.Kp1_b = plastic_modulus1(sp.Ge,sp.R_bar,sp.g_bar,sp.rho1_ratio,h,sp.psi)
            sp.D1 = dilatancy1(sp.R,sp.g,sp.rho1_ratio,sp.psi)
            # print("R g",sp.R,sp.g,sp.rho1_ratio,sp.psi)

        if sp.elastic_flag2 or sp.R == 0.0:
            sp.Kp2_b = 0.0
            sp.D2 = 0.0
        else:
            sign = (sp.p-self.beta)/np.abs(sp.p-self.beta)
            Mg_R = self.M*sp.g/sp.R
            sp.Kp2,sp.Kp2_b = plastic_modulus2(sp.Ge,Mg_R,sp.rho2_ratio,sign)
            sp.D2 = dilatancy2(Mg_R,sign)


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

        if not sp.elastic_flag1:
            theta3_bar = self.Lode_angle(sp.sij)
            dg_bar = self.dg_theta(theta3_bar)
            dF1 = dF1_r(sp.rij_bar,sp.R_bar,theta3_bar,sp.g_bar,dg_bar)
            dF1_tr = np.trace(dF1)
            nij = dF1 - self.I3*dF1_tr/3.0
            sp.nij = nij / np.linalg.norm(nij)
            # print("nij",sp.nij)

        r_abs = math.sqrt(np.square(sp.rij).sum())
        if r_abs == 0.0:
            sp.mij = self.Z3.copy()
        else:
            sp.mij = sp.rij / r_abs         # mij = rij/|rij|
        # print("mij",sp.mij)

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
            # print("Zu",Zu)
            # print("Zd",Zd)
            sp.Zij = Zu / Zd

    # -------------------------------------------------------------------------------------- #
    def set_tensor_Ep(self,sp):
        Lm0 = (np.einsum('pk,ql->pqkl',self.I3,self.I3) + np.einsum('pl,qk->pqkl',self.I3,self.I3))/2.0
        # Lm0 = np.einsum('pk,ql->pqkl',self.I3,self.I3)
        if sp.elastic_flag1:
            Lm1 = np.zeros([3,3,3,3])
        else:
            nD = sp.nij + np.sqrt(2/27)*sp.D1*self.I3
            Lm1 = np.einsum("pq,kl->pqkl",nD,sp.Tij)

        if sp.elastic_flag2:
            Lm2 = np.zeros([3,3,3,3])
        elif sp.R == 0:
            mD = np.sqrt(2/27)*self.I3
            Lm2 = np.einsum("pq,kl->pqkl",mD,sp.Zij)
        else:
            mD = sp.mij + np.sqrt(2/27)*sp.D2*self.I3
            Lm2 = np.einsum("pq,kl->pqkl",mD,sp.Zij)

        Lm = Lm0 - Lm1 - Lm2
        Ee = self.elastic_stiffness(sp.Ge)
        sp.Ep = np.einsum('ijpq,pqkl->ijkl',Ee,Lm)


        # print("Lm",Lm[:,:,0,2])
        # print("Lm0",Lm0[:,:,0,2])
        # print("Lm1",Lm1[:,:,0,2])
        # print("Lm2",Lm2[:,:,0,2])
        #
        # print("Ee_02",Ee[0,2,:,:])
        # print("Ep_02",sp.Ep[:,:,0,2])

        # dep1 = np.einsum('ijkl,kl->ij',Lm1,sp.dstrain)
        # dep2 = np.einsum('ijkl,kl->ij',Lm2,sp.dstrain)
        # print("dstrain",sp.dstrain)
        # print("dep1",dep1)
        # print("dep2",dep2)


    # -------------------------------------------------------------------------------------- #
    def check_unload(self,sp):
        self.set_mapping_stress(sp,False,False)
        self.set_parameters(sp)
        self.set_parameter_nm(sp)
        self.set_parameter_TZ(sp)

        # print("check sp0",sp.elastic_flag1,sp.elastic_flag2)

        dL1 = np.einsum("ij,ij",sp.Tij,sp.dstrain)
        if sp.elastic_flag1 or (dL1 < 0.0):
            elastic_flag1 = True
            self.alpha = np.copy(sp.rij)
        else:
            elastic_flag1 = False
            # print("rho1",sp.rho1_ratio)


        dL2 = np.einsum("ij,ij",sp.Zij,sp.dstrain)
        if sp.elastic_flag2 or (dL2 < -1.e-12):
            elastic_flag2 = True
            self.beta = np.copy(sp.p)
        else:
            elastic_flag2 = False
            # print("rho2",sp.rho2_ratio)


        # print("dL1,dL2",dL1,dL2)
        # print("dstrain",sp.dstrain)
        # print("Tij",sp.Tij)
        # print("Zij",sp.Zij)

        # print("D1",sp.D1,sp.D1*dL1)
        # print("D2",sp.D2,sp.D2*dL2)

        return elastic_flag1,elastic_flag2

    # -------------------------------------------------------------------------------------- #
    def update_parameters(self,sp):
        sp.stress += sp.dstress
        H1,H2 = self.set_mapping_stress(sp)
        self.set_parameters(sp)
        self.set_parameter_nm(sp)
        self.set_parameter_TZ(sp)

        dL1 = np.einsum("ij,ij",sp.Tij,sp.dstrain)
        dL2 = np.einsum("ij,ij",sp.Zij,sp.dstrain)
        # print(" dL:",sp.elastic_flag1,dL1,sp.elastic_flag2,dL2)

        if not sp.elastic_flag1:
            self.L1 += dL1
            self.H1 = H1

        if not sp.elastic_flag2:
            self.H2 = H2


    # -------------------------------------------------------------------------------------- #
    def plastic_stiffness(self,sp):
        self.set_mapping_stress(sp)
        self.set_parameters(sp)
        self.set_parameter_nm(sp)
        self.set_parameter_TZ(sp)
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

    def solve_strain_with_constrain_CU(self,strain_given,stress_given,E,deformation):
        # deformation: True => deform (stress given), False => constrain (strain given)
        n = self.e/(1+self.e)
        Kw = 2.19e9
        Eeps = Kw / n*self.Dijkl

        d = deformation.flatten()
        A = np.reshape(E+Eeps,(9,9))

        strain = np.copy(strain_given.flatten())    #制御するひずみ
        strain[d] = 0.0                        # [0.0,0.0,given,...]deformation Falseでひずみ拘束する場合，ひずみ規定で生じる応力を考える(n+1)

        stress_constrain = np.dot(A,strain)
        stress = np.copy(stress_given.flatten()) - stress_constrain     #n+1でgiven stressとなるように応力を与える

        stress_mask = stress[d]     #deformation Trueの応力の規定により生じる歪を考える

        A_mask = A[d][:,d]

        Ainv = np.linalg.pinv(A_mask)
        strain_mask = Ainv @ stress_mask    #given stressによる変形量

        strain[d] = strain_mask     #constrain_stressによる変形 + 応力規定による変形
        stress = np.dot(A,strain)

        depstress = Kw/n*(strain[0]+strain[4]+strain[8])
        stress -= depstress*np.array([1,0,0,0,1,0,0,0,1])

        return np.reshape(strain,(3,3)), np.reshape(stress,(3,3)), depstress

    # -------------------------------------------------------------------------------------- #
    def elastic_deformation(self,dstrain,dstress,Ge,deformation):
        Ee = self.elastic_stiffness(Ge)
        dstrain_elastic,dstress_elastic = self.solve_strain_with_constrain(dstrain,dstress,Ee,deformation)
        return dstrain_elastic,dstress_elastic

    def plastic_deformation(self,dstrain_given,dstress_given,deformation,sp0,pore_dpressure=0.0,CU=False):
        # sp0 :            stress state defined by EFFECTIVE stress
        # dstress_given :  TOTAL stress

        # extract stress state from sp0
        ef1,ef2 = self.check_unload(sp0)
        strain,eff_stress = sp0.strain,sp0.stress
        pore_pressure     = sp0.pore_pressure
        eff_p = np.trace(eff_stress)/3.0

        # set elasto-plastic tensor in dry condition
        eff_dstress_given = dstress_given - pore_dpressure*np.eye(3)
        sp = self.StateParameters(strain,eff_stress,dstrain_given,eff_dstress_given,self.stress_shift,ef1=ef1,ef2=ef2)
        Ep = self.plastic_stiffness(sp)

        # set Ep tensor in CU condition
        if CU:
            n = self.e/(1+self.e)
            Kw = 2.19e9
            Eeps = Kw / n*self.Dijkl
            Ep += Eeps

        # solve strain and stress increment under the given boundary condition (deformation)
        dstrain_ep,dstress_ep = self.solve_strain_with_constrain(dstrain_given,dstress_given,Ep,deformation)
        if CU:
            pore_dpressure = Kw / n*np.trace(dstrain_ep)
        eff_dstress_ep = dstress_ep - pore_dpressure*np.eye(3)

        # check loading or unloading
        sp_check = self.StateParameters(strain,eff_stress,dstrain_ep,eff_dstress_ep,self.stress_shift,ef1=ef1,ef2=ef2)
        ef1_check,ef2_check = self.check_unload(sp_check)

        if (ef1 != ef1_check) or (ef2 != ef2_check):
            eff_dstress_given = dstress_given - pore_dpressure*np.eye(3)
            sp = self.StateParameters(strain,eff_stress,dstrain_given,eff_dstress_given,self.stress_shift,ef1=ef1_check,ef2=ef2_check)
            Ep = self.plastic_stiffness(sp)
            if CU:
                Ep += Eeps
            dstrain_ep,dstress_ep = self.solve_strain_with_constrain(dstrain_given,dstress_given,Ep,deformation)
            if CU:
                pore_dpressure = Kw / n*np.trace(dstrain_ep)
            eff_dstress_ep = dstress_ep - pore_dpressure*np.eye(3)

        eff_stress_p = eff_stress + eff_dstress_ep
        eff_p_p = np.trace(eff_stress_p)/3
        if eff_p_p + self.stress_shift < self.pmin:
            eff_p_correction = eff_p_p - (self.pmin-self.stress_shift)
            eff_dstress_ep += -eff_p_correction*np.eye(3)
            pore_dpressure += eff_p_correction

        pore_pressure = sp0.pore_pressure + pore_dpressure

        sp2 = self.StateParameters(strain,eff_stress,dstrain_ep,eff_dstress_ep,self.stress_shift,pore_pressure,ef1=ef1_check,ef2=ef2_check)
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
            print("+ N :",ic+1)
            for i in range(nstep):
                print("")
                print("+step:",i,"/",ic)

                dtau = cycle_load(nstep,sr*p0,i)
                dstress_vec = np.array([0.0,0.0,0.0,0.0,0.0,dtau])
                dstress_input = self.vector_to_matrix(dstress_vec)

                sp = self.StateParameters(self.strain,self.stress,sp0.dstrain,dstress_input,self.stress_shift)

                p,R = self.set_stress_variable(self.stress)
                dstrain,dstress,sp0 = \
                    self.plastic_deformation(dstrain_input,dstress_input,deformation,sp)

                # print("stress",self.stress)
                # print("dstress",dstress)
                # print("dstrain",dstrain)

                self.stress += dstress
                self.strain += dstrain

                ev,gamma = self.set_strain_variable(self.strain)
                self.e = self.e0 - ev*(1+self.e0)

                gamma_list += [self.strain[0,2]]
                tau_list += [self.stress[0,2]]
                p_list += [p]
                step_list += [ic*nstep+i]
                ep_list += [(p0-p)/p0]
                # print(p)

                # if p < -3900.0:
                # if ic==3 :
                # if ic==4 and i==3:
                #     plt.plot(p_list,tau_list)
                #     plt.show()
                #     plt.plot(gamma_list,tau_list)
                #     plt.show()
                #     return


        if plot:
            plt.figure()
            plt.plot(gamma_list,tau_list)
            plt.show()
            plt.plot(p_list,tau_list)
            plt.show()
            # plt.plot(step_list,ep_list)
            # plt.show()


    # -------------------------------------------------------------------------------------- #
    def cyclic_shear_test_CU(self,e0,compression_stress,sr=0.2,cycle=20,print_result=False,plot=False):
        def cycle_load(nstep,amp,i):
            tau = amp*np.sin((i+1)/nstep*2*np.pi)
            tau0 = amp*np.sin((i+0)/nstep*2*np.pi)
            return tau - tau0

        self.isotropic_compression(e0,compression_stress)
        # self.stress = np.array([[70e3,0,0],[0,70e3,0],[0,0,70e3]])
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
#
        deformation_vec = np.array([True,True,True,True,True,True],dtype=bool)
        deformation = self.vector_to_matrix(deformation_vec)

        sp0 = self.StateParameters(self.strain,self.stress,dstrain_input,dstress_input,self.stress_shift)

        gamma_list,ev_list,tau_list,ep_list = [],[],[],[]
        p_list,q_list = [],[]
        strain_d, DA_list = [],[]
        epstress_list, epstress_ratio = [],[]
        stressxx,stressyy,stresszz,stresszz_all = [],[],[],[]
        fL_list,h_list,psi_list = [],[],[]
        h1_h2e_list, e_list = [],[]
        H1_list, H2_list, L1_list = [],[],[]

        flag1,flag2,flag5,flag10 = True,True,True,True
        print("sigma_d=",2*sr*p0/1000,"kPa")

        # init effective stress
        self.eff_stress = np.copy(self.stress)
        self.pore_pressure = 0.0
        self.pore_dpressure = 0.0

        for ic in range(ncycle):
            print("###--------------------###")
            print("+ N :",ic+1)
            for i in range(nstep):
                # print("")
                # print("+step:",i,"/",ic)

                dtau = cycle_load(nstep,2*sr*p0,i)
                # print("eps=",eps,"n=",n,"d_ev=",d_ev,"ev",ev)
                dstress_vec = np.array([0.0,0.0,dtau,0.0,0.0,0.0])
                # dstress_vec = np.array([-eps,-eps,dtau-eps,0.0,0.0,0.0])
                dstress_input = self.vector_to_matrix(dstress_vec)

                # set effective stress
                pore_dpressure = self.pore_dpressure
                eff_dstress_input = dstress_input - pore_dpressure*np.eye(3)

                # sp = self.StateParameters(self.strain,self.stress,sp0.dstrain,dstress_input,self.stress_shift)
                sp = self.StateParameters(self.strain,self.eff_stress,sp0.dstrain,eff_dstress_input,self.stress_shift,self.pore_pressure)

                # p,R = self.set_stress_variable(self.stress)
                eff_p,R = self.set_stress_variable(self.eff_stress)

                # dstress: total stress
                dstrain,dstress,sp0 = \
                    self.plastic_deformation(dstrain_input,dstress_input,deformation,sp,pore_dpressure,CU=True)

                self.stress += dstress
                self.strain += dstrain

                self.pore_dpressure = sp0.pore_pressure - self.pore_pressure
                self.pore_pressure = sp0.pore_pressure
                self.eff_stress = self.stress - self.pore_pressure*np.eye(3)
                eff_p = np.trace(self.eff_stress)/3.0
                p = np.trace(self.stress)/3.0

                # print("")
                # print(eff_p,self.eff_stress[2,2]-self.eff_stress[0,0],self.pore_pressure)
                # print(" stress:",self.stress.diagonal())
                # print(" eff_stress:",self.eff_stress.diagonal())

                # if i%50==0:
                #     print("stress",self.stress)
                    # print("strain",self.strain)
                    # print("dstress",dstress)
                    # print("strain",dstrain)
                    # print("fL",self.fL)

                ev,gamma = self.set_strain_variable(self.strain)
                self.e = self.e0 - ev*(1+self.e0)
                dev_strain = self.strain - ev/3.0 * self.I3

                ###---for plot---###
                ep_list += [self.strain[2,2]]
                gamma_list += [gamma]
                tau_list += [self.stress[0,2]]
                ev_list += [ev] #体積ひずみ
                p_list += [eff_p]  #平均有効拘束圧
                q_list += [(self.eff_stress[2,2]-self.eff_stress[0,0])]     #軸差応力
                stressxx += [self.eff_stress[0,0]]
                stressyy += [self.eff_stress[1,1]]
                stresszz += [self.eff_stress[2,2]]
                stresszz_all += [self.stress[2,2]]    #全応力
                epstress_list += [self.pore_pressure]   #過剰間隙水圧
                epstress_ratio += [self.pore_pressure/p]
                strain_d += [self.strain[2,2]]
                fL_list += [self.fL]
                h_list += [self.h]
                psi_list += [self.psi]
                h1_h2e_list += [self.h1-self.h2*self.e]
                e_list += [self.e]
                H1_list += [self.H1]
                H2_list += [self.H2]
                L1_list += [self.L1]

                # if i==0:
                #     exit()
            DA = max(strain_d[-nstep:-1]) - min(strain_d[-nstep:-1])    #各サイクルでの軸ひずみ両振幅
            print("DA",DA*100)
            DA_list += [DA]
            if (1 <= DA*100 < 2) and flag1==True:
                print("1% liquified!!!",max(strain_d[-nstep:-1]),min(strain_d[-nstep:-1]))
                flag1=False
            elif (2 <= DA*100 < 5) and flag2==True:
                print("2% liquified!!!",max(strain_d[-nstep:-1]),min(strain_d[-nstep:-1]))
                flag2=False
            elif (5 <= DA*100 < 10) and flag5==True:
                print("5% liquified!!!",max(strain_d[-nstep:-1]),min(strain_d[-nstep:-1]))
                flag5=False
            elif (10 <= DA*100) and flag10==True:
                print("10% liquified!!!",max(strain_d[-nstep:-1]),min(strain_d[-nstep:-1]))
                flag10=False
            # print("strain",self.strain)
            # print("stress",self.stress)
        step_list = np.linspace(0,cycle,nstep*ncycle)

        gamma_list = np.array(gamma_list)
        tau_list = np.array(tau_list)
        ev_list = np.array(ev_list)
        p_list = np.array(p_list)
        q_list = np.array(q_list)
        stressxx = np.array(stressxx)
        stressyy = np.array(stressyy)
        stresszz = np.array(stresszz)
        stresszz_all = np.array(stresszz_all)
        epstress_list = np.array(epstress_list)
        epstress_ratio = np.array(epstress_ratio)
        strain_d = np.array(strain_d)
        DA_list = np.array(DA_list)
        fL_list =  np.array(fL_list)
        h_list =  np.array(h_list)
        psi_list =  np.array(psi_list)
        h1_h2e_list =  np.array(h1_h2e_list)
        e_list =  np.array(e_list)
        H1_list =  np.array(H1_list)
        H2_list =  np.array(H2_list)
        L1_list =  np.array(L1_list)


        np.savetxt("./result/step.dat",step_list)
        np.savetxt("./result/gamma.dat",gamma_list)
        np.savetxt("./result/tau.dat",tau_list)
        np.savetxt("./result/ev.dat",ev_list)
        np.savetxt("./result/p.dat",p_list)
        np.savetxt("./result/q.dat",q_list)
        np.savetxt("./result/stressxx.dat",stressxx)
        np.savetxt("./result/stressyy.dat",stressyy)
        np.savetxt("./result/stresszz.dat",stresszz)
        np.savetxt("./result/stresszz_all.dat",stresszz_all)
        np.savetxt("./result/epstress.dat",epstress_list)
        np.savetxt("./result/epstrassratio.dat",epstress_ratio)
        np.savetxt("./result/strain_d.dat",strain_d)
        np.savetxt("./result/DA.dat",DA_list)
        np.savetxt("./result/fL.dat",fL_list)
        np.savetxt("./result/h.dat",h_list)
        np.savetxt("./result/psi.dat",psi_list)
        np.savetxt("./result/h1-h2e.dat",h1_h2e_list)
        np.savetxt("./result/e.dat",e_list)
        np.savetxt("./result/H1.dat",H1_list)
        np.savetxt("./result/H2.dat",H2_list)
        np.savetxt("./result/L1.dat",L1_list)

        if plot:
            plt.figure()
            plt.plot(gamma_list,tau_list)
            plt.xlabel("gamma")
            plt.ylabel("tau")
            plt.show()
            plt.plot(p_list,tau_list)
            plt.xlabel("p")
            plt.ylabel("tau")
            plt.show()
            plt.plot(gamma_list,ep_list)
            plt.show()





    # -------------------------------------------------------------------------------------- #
    def cyclic_pure_shear_test(self,e0,compression_stress,gmax=0.01,cycle=30,print_result=False,plot=False):
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

        print(self.strain)
        print(self.stress)
        print(self.e)
        self.clear_strain()
        ###

        nstep = 4000
        ncycle = cycle

        dstrain_vec = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        dstrain_input = self.vector_to_matrix(dstrain_vec)

        dstress_vec = np.array([0.0,0.0,0.0,0.0,0.0,0.0])
        dstress_input = self.vector_to_matrix(dstress_vec)

        deformation_vec = np.array([False,False,True,False,False,False],dtype=bool)
        # deformation_vec = np.array([False,False,False,False,False,False],dtype=bool)
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
                ep_list += [self.strain[2,2]]

                # print(self.strain)
                # print(self.stress)

            print(self.stress[0,0],self.stress[1,1],self.stress[2,2])
            print(self.e,self.strain[2,2])

        if plot:
            plt.figure()
            plt.plot(gamma_list,tau_list)
            plt.xlabel("gamma")
            plt.ylabel("tau")
            plt.show()
            plt.plot(p_list,tau_list)
            plt.xlabel("p")
            plt.ylabel("tau")
            plt.show()
            plt.plot(gamma_list,ep_list)
            plt.show()


# --------------------------------#
if __name__ == "__main__":

    # e0 = 0.83
    # Li_model = Li2002(G0=60,nu=0.33,M=1.61,c=0.75,eg=0.99,rlambdac=0.019,xi=0.7, \
    #              d1=0.41,m=3.5,h1=2.1,h2=2.03,h3=2.2,n=1.1, \
    #              d2=1,h4=3.5,a=1,b1=0.005,b2=2.0,b3=0.01,cohesion=0.0,e0=e0)

    e0 = 0.85
    Li_model = Li2002(G0=60,nu=0.33,M=1.61,c=0.75,eg=0.99,rlambdac=0.019,xi=0.7, \
                 d1=0.41,m=3.5,h1=2.1,h2=2.03,h3=2.2,n=1.1, \
                 d2=1,h4=3.5,a=1,b1=0.005,b2=2.0,b3=0.01,cohesion=0.0,e0=e0)


    compression_stress = 70e3
    # Li_model.cyclic_shear_test_CU(e0,compression_stress,sr=0.25,cycle=30,print_result=True,plot=False)
    Li_model.cyclic_shear_test_CU(e0,compression_stress,sr=0.1,cycle=30,print_result=True,plot=False)
    exit()



    cs_list = [1.e3,2.5e3,5.e3,10.e3,20.e3,40.e3,300.e3]
    sigma1_list, sigma3_list = [],[]
    gamma_list, ev_list = [],[]
    for compression_stress in cs_list:
        # Li_model = Li2002(G0=420,nu=0.33,M=0.97,eg=0.957,d1=0.0)
        Li_model = Li2002(G0=210,nu=0.33,M=0.93,eg=0.957,d1=0.41,cohesion=4.e3)
        # Li_model = Li2002()
        s1,s3,gamma,ev = Li_model.triaxial_compression(e0,compression_stress,print_result=True,plot=False)

        sigma1_list += [s1]
        sigma3_list += [s3]
        gamma_list += [gamma]
        ev_list += [ev]

    # -----------------
    phi = 30.0  # [deg]
    c = 3.e3    # [Pa]
    beta = 2.5  # [deg]

    x = np.linspace(0,np.max(sigma1_list),100)
    y = c + x*np.tan(np.deg2rad(phi))


    fig = plt.figure()
    # ax = plt.axes()
    # plt.grid()

    for s1,s3 in zip(sigma1_list,sigma3_list):
        xc = 0.5*(s1+s3)
        r = 0.5*(s1-s3)
        c = patches.Circle(xy=(xc,0),radius=r,fill=False,ec='k')
        # ax.add_patch(c)

    # plt.plot(x,y)
    #
    # plt.axis('scaled')
    # ax.set_aspect('equal')
    #
    # plt.show()

    # -----------------
    x = np.array(gamma_list[0])
    y = -x*np.tan(np.deg2rad(beta))

    plt.grid()
    for gamma,ev in zip(gamma_list,ev_list):
        plt.plot(gamma,ev,c="k")
    # plt.plot(x,y)
    plt.axis('auto')
    plt.show()
