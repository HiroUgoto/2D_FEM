import numpy as np
import torch as t
import os,pickle
import matplotlib.pyplot as plt
from collections import deque

from ep_model.DL1d import numerical_model as nm
from ep_model.DL1d import Net

def load(name):
    dirname = os.path.dirname(__file__)
    fname = os.path.join(dirname,'save',name)
    with open(fname,'rb') as f:
        obj = pickle.load(f)
    return obj


class DL1d:
    def __init__(self,info,savelen=None,maxlen=None):
        '''
        Args:
            `info (dict)`: Dictionary that contains material parameters. Keys are 'H', 'P0', 'N', 'G0', 'sand', 'silt', 'clay', 'wL' and 'wP'. The first four are in SI units and the last five are in percent(%).
        '''
        self.info_dict = info
        if maxlen is None:
            self.enc_gamma = [0]*100
            self.enc_nmtau1 = [0]*100
            self.enc_nmtau2 = [0]*100
        else:
            self.enc_gamma = deque([0]*100,maxlen)
            self.enc_nmtau1 = deque([0]*100,maxlen)
            self.enc_nmtau2 = deque([0]*100,maxlen)
        if savelen is None:
            self.gamma_list = []
            self.tau_list = []
        else:
            self.gamma_list = deque([],savelen)
            self.tau_list = deque([],savelen)
        self.gamma = 0.0
        self.tau = 0.0

        self._init_encode()
        self._init_dl_models()

    def _init_encode(self):
        # info
        self.f_name = ['N','H','G0','sand','silt','clay','wL','wP','P0']
        mean_dict,std_dict = load('info_mean_std.pickle')
        self.imean = np.array([mean_dict[key] for key in self.f_name])
        self.istd = np.array([std_dict[key] for key in self.f_name])

        # strain, stress
        self.refgamma = 0.05/5
        self.reftau = 1e6/5

    def _init_dl_models(self):
        self.dl_gr = Net.DL_gr()
        self.dl_h = Net.DL_h()
        self.dl_s2s = Net.DL_s2s()

    def initial_state(self,info_update):
        self.info_dict.update(info_update)
        info = np.array([self.info_dict[key] for key in self.f_name])
        self.info_enc = self._encode_info(info)
        self._init_soil()

    def _init_soil(self):
        '''
        Args:
            ``info_enc``: torch.Tensor of shape (nfeature,)
        '''
        gr = self.dl_gr(self.info_enc)[0]
        hmax = self.dl_h(self.info_enc)[0]
        # print(f'gr {gr:.2e}, G0 {self.G0:.2e}, hmax {hmax:.2f}')
        self.ghe_param = {
            'gr': gr,
            'G0': self.info_dict['G0'],
            'hmax': hmax,
            'c1_0': 1,
            'c1_inf': 0.170,
            'c2_0': 0.830,
            'c2_inf': 2.5,
            'alpha': 2.860,
            'beta': 3.229,
        }
        self._init_numerical_models()

    def _init_numerical_models(self):
        self.nm1 = nm.mod_GHE(**self.ghe_param)
        self.nm2 = nm.GHE_S(**self.ghe_param)


    def _encode_info(self,info):
        x = (info-self.imean)/self.istd
        x = np.nan_to_num(x,nan=-100)
        return t.from_numpy(x.astype(np.float32))

    def _encode_gamma(self,gamma):
        x = np.tanh(gamma/self.refgamma)
        return x

    def _encode_tau(self,tau):
        x = np.tanh(tau/self.reftau)
        return x

    def _decode_tau(self,tau):
        x = t.atanh(tau)*self.reftau
        return x.detach().numpy()


    def shear(self,gamma):
        '''
        Args:
            `gamma (float)`: input gamma.
        '''
        nmtau1 = self.nm1.shear(gamma)
        nmtau2 = self.nm2.shear(gamma)

        self.gamma = gamma
        self.enc_gamma.append(self._encode_gamma(gamma))
        self.enc_nmtau1.append(self._encode_tau(nmtau1))
        self.enc_nmtau2.append(self._encode_tau(nmtau2))

        enc_gamma = t.Tensor(self.enc_gamma)
        enc_nmtau = t.Tensor([self.enc_nmtau1,self.enc_nmtau2]).T
        x = self.info_enc,enc_gamma,enc_nmtau
        with t.no_grad():
            enc_tau = t.squeeze(self.dl_s2s(*x))
        self.tau = self._decode_tau(enc_tau[-1])
        self.gamma_list.append(self.gamma)
        self.tau_list.append(self.tau)
        return self.tau

    def shear_d(self,dgamma):
        '''
        Args:
            `dgamma (float)`: gamma[t]-gamma[t-1]
        '''
        gamma = self.gamma+dgamma
        tau0 = self.tau
        tau = self.shear(gamma)
        return tau-tau0,gamma

    def shear_nm(self,gamma):
        '''Returns stress from DL model and numerical models.'''
        nmtau1 = self.nm1.shear(gamma)
        nmtau2 = self.nm2.shear(gamma)

        self.gamma = gamma
        self.enc_gamma.append(self._encode_gamma(gamma))
        self.enc_nmtau1.append(self._encode_tau(nmtau1))
        self.enc_nmtau2.append(self._encode_tau(nmtau2))

        enc_gamma = t.Tensor(self.enc_gamma)
        enc_nmtau = t.Tensor([self.enc_nmtau1,self.enc_nmtau2]).T
        x = self.info_enc,enc_gamma,enc_nmtau
        with t.no_grad():
            enc_tau = t.squeeze(self.dl_s2s(*x))
        self.tau = self._decode_tau(enc_tau[-1])
        self.gamma_list.append(self.gamma)
        self.tau_list.append(self.tau)
        return self.tau,nmtau1,nmtau2

    def cyclic_shear_test(self,gamma=None,maxgamma=0.02,plot=False):
        if gamma is None:
            cycle = 10
            steppercycle = 300
            step0,step1 = steppercycle,cycle*steppercycle
            tim1 = np.linspace(0,cycle,step1)
            mask1 = np.ones_like(tim1)
            n1,n2 = int(step1/3),int(step1*2/3)
            mask1[:n1] = np.linspace(0,1,n1)
            mask1[n2:] = np.linspace(1,0,step1-n2)
            gamma1 = np.sin(tim1*2*np.pi) * mask1 * maxgamma
            gamma = np.append(np.zeros(step0),gamma1)

        nseq = len(gamma)
        tau_list = []
        nmtau1_list = []
        nmtau2_list = []
        for i in range(nseq):
            tau,nmtau1,nmtau2 = self.shear_nm(gamma[i])
            tau_list.append(tau)
            nmtau1_list.append(nmtau1)
            nmtau2_list.append(nmtau2)

        if plot:
            plt.figure()
            plt.plot(gamma,nmtau1_list,label='nm1')
            plt.plot(gamma,nmtau2_list,label='nm2')
            plt.plot(gamma,tau_list,label='dl')
            plt.legend()
            plt.show()

        return tau_list

    def plot(self,fname='result/constitution.png'):
        gamma = 100*np.array(self.gamma_list)
        tau = 1/1000*np.array(self.tau_list)
        fig,ax = plt.subplots(figsize=[5,5])
        ax.set_xlabel('gamma [%]')
        ax.set_ylabel('tau [kN/m2]', labelpad=4.0)
        ax.plot(gamma,tau) #label
        fig.savefig(fname)
        plt.close(fig)


class Light_DL1d(DL1d):
    def __init__(self, info, maxlen=500):
        '''
        Args:
            `info (dict)`: Dictionary that contains material parameters. Keys are 'H', 'P0', 'N', 'G0', 'sand', 'silt', 'clay', 'wL' and 'wP'. The first four are in SI units and the last five are in percent(%).
        '''
        self.G0 = info['G0']
        self.enc_gamma = deque([0]*100,maxlen)
        self.enc_nmtau1 = deque([0]*100,maxlen)
        self.enc_nmtau2 = deque([0]*100,maxlen)
        self.gamma_list = []
        self.tau_list = []
        self.gamma = 0.0
        self.tau = 0.0

        self._init_encode(info)
        self._init_dl_models()
        self._init_soil()
        self._init_numerical_models()


class GHE_mix(DL1d):
    def __init__(self,info):
        self.G0 = info['G0']
        self.gamma_list = []
        self.tau_list = []
        self.gamma = 0.0
        self.tau = 0.0

        print(11,info)
        self._init_encode(info)
        self._init_dl_models()
        self._init_soil()
        self._init_numerical_models()

    def _init_dl_models(self):
        self.dl_gr = Net.DL_gr()
        self.dl_h = Net.DL_h()

    def shear(self,gamma,rate=0.5):
        tau1 = self.nm1.shear(gamma)
        tau2 = self.nm2.shear(gamma)

        self.gamma = gamma
        self.tau = rate*tau1 + (1-rate)*tau2
        self.gamma_list.append(self.gamma)
        self.tau_list.append(self.tau)
        return self.tau

    def shear_nm(self,gamma,rate=0.5):
        tau1 = self.nm1.shear(gamma)
        tau2 = self.nm2.shear(gamma)

        self.gamma = gamma
        self.tau = rate*tau1 + (1-rate)*tau2
        self.gamma_list.append(self.gamma)
        self.tau_list.append(self.tau)
        return self.tau,tau1,tau2

    def plot(self,fname='result/constitution_ghe.png'):
        pass


if __name__ == "__main__":
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
    model = DL1d(info)
    model.cyclic_shear_test()
