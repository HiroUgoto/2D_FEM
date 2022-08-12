import numpy as np

class Source:
    def __init__(self,dip,width,id):
        self.element_id = id
        self.dip = np.deg2rad(dip)
        self.width = width

        self.set_strain_tensor()


    def set_strain_tensor(self):
        self.strain_tensor = np.zeros(3,dtype=np.float64)

        self.strain_tensor[0] = -np.sin(2*self.dip)
        self.strain_tensor[1] =  np.sin(2*self.dip)
        self.strain_tensor[2] =  np.cos(2*self.dip)
