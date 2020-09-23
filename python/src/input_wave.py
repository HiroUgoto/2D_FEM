import numpy as np

def ricker(tim,tp,fp,amp):
    t1 = ((tim-tp)*np.pi*fp)**2
    return (2*t1-1)*np.exp(-t1)*amp
