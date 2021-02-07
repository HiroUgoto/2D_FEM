import numpy as np

class Node:
    def __init__(self,id,xyz,freedom):
        self.id = id
        self.xyz = xyz
        self.freedom = freedom
        self.dof = len(freedom)

    def print(self):
        print(self.id,":",self.xyz,",",self.freedom)

    def set_initial_condition(self):
        self.u  = np.zeros(self.dof,dtype=np.float64)
        self.um = np.zeros(self.dof,dtype=np.float64)
        self.v  = np.zeros(self.dof,dtype=np.float64)

        self.mass   = np.zeros(self.dof,dtype=np.float64)
        self.c      = np.zeros(self.dof,dtype=np.float64)
        self.k      = np.zeros(self.dof,dtype=np.float64)

        self.force = np.zeros(self.dof,dtype=np.float64)
        self.static_force = np.zeros(self.dof,dtype=np.float64)
        self.dynamic_force = np.zeros(self.dof,dtype=np.float64)

        self._up = np.zeros(self.dof,dtype=np.float64)
        self._ur = np.zeros(self.dof,dtype=np.float64)
        self._uy = np.zeros(self.dof,dtype=np.float64)
