import numpy as np

def set_source(dof,elements,dip,width,sx,sz,n=1):
    source_list = []
    dip_rad = np.deg2rad(dip)

    wn = np.linspace(-width/2,width/2,n+1)
    w = np.convolve(wn,[0.5,0.5],"valid")
    dw = width/n

    x = np.array([sx,sz])
    id = 0
    for i in range(n):
        x[0] = sx + w[i]*np.cos(dip_rad)
        x[1] = sz + w[i]*np.sin(dip_rad)

        for element in elements:
            if element.dim == 2:
                is_inside,xi = element.check_inside(x)
                if is_inside:
                    source = Source(dof,i,dip,dw,element.id,xi[0],xi[1])
                    source_list += [source]
                    id += 1
                    break

    return source_list


class Source:
    def __init__(self,dof,id,dip,width,element_id,xi,zeta):
        self.dof = dof
        self.id = id
        self.element_id = element_id
        self.xi,self.zeta = xi, zeta

        self.dip = np.deg2rad(dip)
        self.width = width

        self.set_strain_tensor()

    def print(self):
        print(self.id,":",self.dip,",",self.width)
        print("    ",self.element_id,",",(self.xi,self.zeta))

    def set_strain_tensor(self):
        if self.dof == 1:
            self.strain_tensor = np.zeros(2,dtype=np.float64)

            self.strain_tensor[0] =  np.sin(self.dip) *self.width
            self.strain_tensor[1] = -np.cos(self.dip) *self.width

        elif self.dof == 2:
            self.strain_tensor = np.zeros(3,dtype=np.float64)

            self.strain_tensor[0] = -np.sin(2*self.dip) *self.width
            self.strain_tensor[1] =  np.sin(2*self.dip) *self.width
            self.strain_tensor[2] =  np.cos(2*self.dip) *self.width
