from functools import lru_cache
import numpy as np

def set_style(style):
    if style == "2d4solid":
        return Solid_2d_4Node()
    elif style == "2d8solid":
        return Solid_2d_8Node()
    elif style == "2d9solid":
        return Solid_2d_9Node()
    elif style == "1d2line":
        return Line_1d_2Node()
    elif style == "1d3line":
        return Line_1d_3Node()
    elif style == "1d2input":
        return Input_1d_2Node()
    elif style == "1d3input":
        return Input_1d_3Node()
    elif style == "connect":
        return Connect()

# =================== Element style classes ============================ #
class Gauss_Points:
    def __init__(self,dn,w,N,M=None):
        self.dn = dn
        self.w = w
        self.N = N
        self.M = M

# =================== Element style classes ============================ #
class Solid_2d_4Node:
    def __init__(self):
        self.dim = 2
        self.gauss = np.polynomial.legendre.leggauss(3)

    def init_dn(self,n):
        return np.zeros([n,n,4,2])

    @lru_cache()
    def shape_function_n(self,xi,zeta):
        n = np.zeros(4)
        n[0] = (1.0 - xi)*(1.0 - zeta) / 4.0
        n[1] = (1.0 + xi)*(1.0 - zeta) / 4.0
        n[2] = (1.0 + xi)*(1.0 + zeta) / 4.0
        n[3] = (1.0 - xi)*(1.0 + zeta) / 4.0
        return n

    @lru_cache()
    def shape_function_dn(style,xi,zeta):
        dn = np.zeros([4,2])
        dn[0,0] = -(1.0 - zeta) / 4.0
        dn[0,1] = -(1.0 -   xi) / 4.0

        dn[1,0] =  (1.0 - zeta) / 4.0
        dn[1,1] = -(1.0 +   xi) / 4.0

        dn[2,0] =  (1.0 + zeta) / 4.0
        dn[2,1] =  (1.0 +   xi) / 4.0

        dn[3,0] = -(1.0 + zeta) / 4.0
        dn[3,1] =  (1.0 -   xi) / 4.0
        return dn

# ---------------------------------------------------------------------- #
class Solid_2d_8Node:
    def __init__(self):
        self.dim = 2
        self.gauss = np.polynomial.legendre.leggauss(5)

    def init_dn(self,n):
        return np.zeros([n,n,8,2])

    @lru_cache()
    def shape_function_n(self,xi,zeta):
        n = np.zeros(8)
        n[0] = (1.0 - xi)*(1.0 - zeta)*(-1.0-xi-zeta) / 4.0
        n[1] = (1.0 + xi)*(1.0 - zeta)*(-1.0+xi-zeta) / 4.0
        n[2] = (1.0 + xi)*(1.0 + zeta)*(-1.0+xi+zeta) / 4.0
        n[3] = (1.0 - xi)*(1.0 + zeta)*(-1.0-xi+zeta) / 4.0

        n[4] = (1.0 - xi*xi)*(1.0 - zeta) / 2.0
        n[5] = (1.0 + xi)*(1.0 - zeta*zeta) / 2.0
        n[6] = (1.0 - xi*xi)*(1.0 + zeta) / 2.0
        n[7] = (1.0 - xi)*(1.0 - zeta*zeta) / 2.0
        return n

    @lru_cache()
    def shape_function_dn(style,xi,zeta):
        dn = np.zeros([8,2])
        dn[0,0] = (1.0 - zeta)*(2.0*xi+zeta) / 4.0
        dn[0,1] = (1.0 -   xi)*(xi+2.0*zeta) / 4.0

        dn[1,0] = (1.0 - zeta)*(2.0*xi-zeta) / 4.0
        dn[1,1] = -(1.0 +  xi)*(xi-2.0*zeta) / 4.0

        dn[2,0] = (1.0 + zeta)*(2.0*xi+zeta) / 4.0
        dn[2,1] = (1.0 +   xi)*(xi+2.0*zeta) / 4.0

        dn[3,0] = (1.0 + zeta)*(2.0*xi-zeta) / 4.0
        dn[3,1] = -(1.0 -  xi)*(xi-2.0*zeta) / 4.0

        dn[4,0] = -xi*(1.0 - zeta)
        dn[4,1] = (xi**2-1.0) / 2.0

        dn[5,0] = (1.0 - zeta**2) / 2.0
        dn[5,1] = -(1.0 + xi)*zeta

        dn[6,0] = -xi*(1.0 + zeta)
        dn[6,1] = (1.0 - xi**2) / 2.0

        dn[7,0] = -(1.0 - zeta**2) / 2.0
        dn[7,1] = -(1.0 - xi)*zeta
        return dn

# ---------------------------------------------------------------------- #
class Solid_2d_9Node:
    def __init__(self):
        self.dim = 2
        self.gauss = np.polynomial.legendre.leggauss(5)

    def init_dn(self,n):
        return np.zeros([n,n,9,2])

    @lru_cache()
    def shape_function_n(self,xi,zeta):
        n = np.zeros(9)
        n[0] =  (1.0 - xi)*(1.0 - zeta)*xi*zeta / 4.0
        n[1] = -(1.0 + xi)*(1.0 - zeta)*xi*zeta / 4.0
        n[2] =  (1.0 + xi)*(1.0 + zeta)*xi*zeta / 4.0
        n[3] = -(1.0 - xi)*(1.0 + zeta)*xi*zeta / 4.0

        n[4] = -(1.0 - xi*xi)*zeta*(1.0 - zeta) / 2.0
        n[5] =  (1.0 + xi)*xi*(1.0 - zeta*zeta) / 2.0
        n[6] =  (1.0 - xi*xi)*zeta*(1.0 + zeta) / 2.0
        n[7] = -(1.0 - xi)*xi*(1.0 - zeta*zeta) / 2.0

        n[8] = (1.0-xi*xi)*(1.0-zeta*zeta)
        return n

    @lru_cache()
    def shape_function_dn(style,xi,zeta):
        dn = np.zeros([9,2])
        dn[0,0] =  (2.0*xi-1.0)*(zeta-1.0)*zeta / 4.0
        dn[0,1] =  (xi-1.0)*xi*(2.0*zeta-1.0) / 4.0

        dn[1,0] =  (2.0*xi+1.0)*(zeta-1.0)*zeta / 4.0
        dn[1,1] =  (xi+1.0)*xi*(2.0*zeta-1.0) / 4.0

        dn[2,0] =  (2.0*xi+1.0)*(zeta+1.0)*zeta / 4.0
        dn[2,1] =  (xi+1.0)*xi*(2.0*zeta+1.0) / 4.0

        dn[3,0] =  (2.0*xi-1.0)*(zeta+1.0)*zeta / 4.0
        dn[3,1] =  (xi-1.0)*xi*(2.0*zeta+1.0) / 4.0

        dn[4,0] =  xi*(1.0-zeta)*zeta
        dn[4,1] =  (1.0-xi*xi)*(2.0*zeta-1.0) / 2.0

        dn[5,0] =  (2.0*xi+1.0)*(1.0-zeta*zeta) / 2.0
        dn[5,1] =  -xi*(1.0+xi)*zeta

        dn[6,0] =  -xi*(1.0+zeta)*zeta
        dn[6,1] =  (1.0-xi*xi)*(2.0*zeta+1.0) / 2.0

        dn[7,0] =  (2.0*xi-1.0)*(1.0-zeta*zeta) / 2.0
        dn[7,1] =  xi*(1.0-xi)*zeta

        dn[8,0] =  -2.0*xi*(1.0-zeta*zeta)
        dn[8,1] =  -2.0*(1.0-xi*xi)*zeta
        return dn

# ---------------------------------------------------------------------- #
class Line_1d_2Node:
    def __init__(self):
        self.dim = 1
        self.gauss = np.polynomial.legendre.leggauss(3)

    def init_dn(self,n):
        return np.zeros([n,2])

    @lru_cache()
    def shape_function_n(self,xi,zeta=0.0):
        n = np.zeros(2)
        n[0] = (1.0 - xi) / 2.0
        n[1] = (1.0 + xi) / 2.0
        return n

    @lru_cache()
    def shape_function_dn(style,xi,zeta=0.0):
        dn = np.zeros(2)
        dn[0] = -0.5
        dn[1] =  0.5
        return dn

# ---------------------------------------------------------------------- #
class Line_1d_3Node:
    def __init__(self):
        self.dim = 1
        self.gauss = np.polynomial.legendre.leggauss(5)

    def init_dn(self,n):
        return np.zeros([n,3])

    @lru_cache()
    def shape_function_n(self,xi,zeta=0.0):
        n = np.zeros(3)
        n[0] = -xi*(1.0 - xi) / 2.0
        n[1] =  xi*(1.0 + xi) / 2.0
        n[2] = (1.0 - xi)*(1.0 + xi)
        return n

    @lru_cache()
    def shape_function_dn(style,xi,zeta=0.0):
        dn = np.zeros(3)
        dn[0] = xi - 0.5
        dn[1] = xi + 0.5
        dn[2] = -2.0*xi
        return dn

# ---------------------------------------------------------------------- #
class Connect:
    def __init__(self):
        self.dim = 0
        self.gauss = np.array([]),np.array([])

# ---------------------------------------------------------------------- #

class Input_1d_2Node(Line_1d_2Node):
    pass

# ---------------------------------------------------------------------- #
class Input_1d_3Node(Line_1d_3Node):
    pass
