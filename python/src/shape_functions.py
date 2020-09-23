import numpy as np

def shape_function_n(style,xi,zeta):
    n = np.zeros(style)

    if style == 4:
        n[0] = (1.0 - xi)*(1.0 - zeta) / 4.0
        n[1] = (1.0 + xi)*(1.0 - zeta) / 4.0
        n[2] = (1.0 + xi)*(1.0 + zeta) / 4.0
        n[3] = (1.0 - xi)*(1.0 + zeta) / 4.0

    elif style == 8:
        n[0] = (1.0 - xi)*(1.0 - zeta)*(-1.0-xi-zeta) / 4.0
        n[1] = (1.0 + xi)*(1.0 - zeta)*(-1.0+xi-zeta) / 4.0
        n[2] = (1.0 + xi)*(1.0 + zeta)*(-1.0+xi+zeta) / 4.0
        n[3] = (1.0 - xi)*(1.0 + zeta)*(-1.0-xi+zeta) / 4.0

        n[4] = (1.0 - xi*xi)*(1.0 - zeta) / 2.0
        n[5] = (1.0 + xi)*(1.0 - zeta*zeta) / 2.0
        n[6] = (1.0 - xi*xi)*(1.0 + zeta) / 2.0
        n[7] = (1.0 - xi)*(1.0 - zeta*zeta) / 2.0

    elif style == 9:
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

def shape_function_dn(style,xi,zeta):
    dn = np.zeros([style,2])

    if style == 4:
        dn[0,0] = -(1.0 - zeta) / 4.0
        dn[0,1] = -(1.0 -   xi) / 4.0

        dn[1,0] =  (1.0 - zeta) / 4.0
        dn[1,1] = -(1.0 +   xi) / 4.0

        dn[2,0] =  (1.0 + zeta) / 4.0
        dn[2,1] =  (1.0 +   xi) / 4.0

        dn[3,0] = -(1.0 + zeta) / 4.0
        dn[3,1] =  (1.0 -   xi) / 4.0

    elif style == 8:
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

    elif style == 9:
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
