import numpy as np

def set_dim(style):
    if "1d" in style:
        return 1
    elif "2d" in style:
        return 2

def set_gauss(style):
    if style == "2d4solid":
        return np.polynomial.legendre.leggauss(3)
    elif style == "2d8solid":
        return np.polynomial.legendre.leggauss(5)
    elif style == "2d9solid":
        return np.polynomial.legendre.leggauss(5)
    elif style == "1d2line":
        return np.polynomial.legendre.leggauss(3)
    elif style == "1d3line":
        return np.polynomial.legendre.leggauss(5)
    elif style == "1d2input":
        return np.polynomial.legendre.leggauss(3)
    elif style == "1d3input":
        return np.polynomial.legendre.leggauss(5)

def shape_function_n(style,xi,zeta=0.0):
    if style == "1d2line":
        n = np.zeros(2)
        n[0] = (1.0 - xi) / 2.0
        n[1] = (1.0 + xi) / 2.0

    elif style == "1d3line":
        n = np.zeros(3)
        n[0] = -xi*(1.0 - xi) / 2.0
        n[1] =  xi*(1.0 + xi) / 2.0
        n[2] = (1.0 - xi)*(1.0 + xi)

    elif style == "1d2input":
        n = np.zeros(2)
        n[0] = (1.0 - xi) / 2.0
        n[1] = (1.0 + xi) / 2.0

    elif style == "1d3input":
        n = np.zeros(3)
        n[0] = -xi*(1.0 - xi) / 2.0
        n[1] =  xi*(1.0 + xi) / 2.0
        n[2] = (1.0 - xi)*(1.0 + xi)

    elif style == "2d4solid":
        n = np.zeros(4)
        n[0] = (1.0 - xi)*(1.0 - zeta) / 4.0
        n[1] = (1.0 + xi)*(1.0 - zeta) / 4.0
        n[2] = (1.0 + xi)*(1.0 + zeta) / 4.0
        n[3] = (1.0 - xi)*(1.0 + zeta) / 4.0

    elif style == "2d8solid":
        n = np.zeros(8)
        n[0] = (1.0 - xi)*(1.0 - zeta)*(-1.0-xi-zeta) / 4.0
        n[1] = (1.0 + xi)*(1.0 - zeta)*(-1.0+xi-zeta) / 4.0
        n[2] = (1.0 + xi)*(1.0 + zeta)*(-1.0+xi+zeta) / 4.0
        n[3] = (1.0 - xi)*(1.0 + zeta)*(-1.0-xi+zeta) / 4.0

        n[4] = (1.0 - xi*xi)*(1.0 - zeta) / 2.0
        n[5] = (1.0 + xi)*(1.0 - zeta*zeta) / 2.0
        n[6] = (1.0 - xi*xi)*(1.0 + zeta) / 2.0
        n[7] = (1.0 - xi)*(1.0 - zeta*zeta) / 2.0

    elif style == "2d9solid":
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

def shape_function_dn(style,xi,zeta=0.0):
    if style == "1d2line":
        dn = np.zeros(2)
        dn[0] = -0.5
        dn[1] =  0.5

    elif style == "1d3line":
        dn = np.zeros(3)
        dn[0] = xi - 0.5
        dn[1] = xi + 0.5
        dn[2] = -2.0*xi

    elif style == "1d2input":
        dn = np.zeros(2)
        dn[0] = -0.5
        dn[1] =  0.5

    elif style == "1d3input":
        dn = np.zeros(3)
        dn[0] = xi - 0.5
        dn[1] = xi + 0.5
        dn[2] = -2.0*xi

    elif style == "2d4solid":
        dn = np.zeros([4,2])
        dn[0,0] = -(1.0 - zeta) / 4.0
        dn[0,1] = -(1.0 -   xi) / 4.0

        dn[1,0] =  (1.0 - zeta) / 4.0
        dn[1,1] = -(1.0 +   xi) / 4.0

        dn[2,0] =  (1.0 + zeta) / 4.0
        dn[2,1] =  (1.0 +   xi) / 4.0

        dn[3,0] = -(1.0 + zeta) / 4.0
        dn[3,1] =  (1.0 -   xi) / 4.0

    elif style == "2d8solid":
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

    elif style == "2d9solid":
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
