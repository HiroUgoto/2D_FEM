import numpy as np

def ricker(tim,fp,tp,amp):
    t1 = ((tim-tp)*np.pi*fp)**2
    return (2*t1-1)*np.exp(-t1)*amp


def tapered_sin(tim,fp,taper,duration,amp):
    wave = simple_sin(tim,fp,amp)
    coeff = np.ones_like(wave)

    ind = np.where((0 <= tim) & (tim < taper))
    coeff[ind] = tim[ind]/taper

    ind = np.where((duration-taper < tim) & (tim <= duration))
    coeff[ind] = (duration-tim[ind])/taper

    ind = np.where(duration < tim)
    coeff[ind] = 0.0

    return wave*coeff


def simple_sin(tim,fp,amp):
    return amp*np.sin(2*np.pi*fp*tim)


def smoothed_ramp(tim,fp,tp,amp=1.0):
    return amp * (1.+np.tanh(4*fp*(tim-tp)))/2.
