import numpy as np
from ep_model.DL1d import constitution

# Using soil 10 of tabdata for an example.
info = {
    'H':16.0,
    'P0':108,
    'N':4.0,
    'G0':26.2e6,
    'sand':4.0,
    'silt':40.0,
    'clay':56.0,
    'wL':np.nan,
    'wP':np.nan,
}
model = constitution.DL1d(info)
model.cyclic_shear_test()