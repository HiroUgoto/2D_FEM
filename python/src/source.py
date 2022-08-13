import numpy as np

def set_source(elements,dip,width,sx,sz,n=1):
    source_list = []
    dip_rad = np.deg2rad(dip)

    wn = np.linspace(-width/2,width/2,n+1)
    w = np.convolve(wn,[0.5,0.5],"valid")
    dw = width/n

    x = np.array([sx,sz])
    for i in range(n):
        x[0] = sx + w[i]*np.cos(dip_rad)
        x[1] = sz + w[i]*np.sin(dip_rad)

        for element in elements:
            if element.dim == 2:
                is_inside,xi = element.check_inside(x)
                if is_inside:
                    source = Source(i,dip,dw,element.id,xi[0],xi[1])
                    source_list += [source]

    return source_list


class Source:
    def __init__(self,id,dip,width,element_id,xi,zeta):
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
        self.strain_tensor = np.zeros(3,dtype=np.float64)

        self.strain_tensor[0] = -np.sin(2*self.dip) *self.width
        self.strain_tensor[1] =  np.sin(2*self.dip) *self.width
        self.strain_tensor[2] =  np.cos(2*self.dip) *self.width
