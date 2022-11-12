import os
import numpy as np
import torch as t
import torch.nn as nn
from ep_model.DL1d import TCN_class as TCN

dirname = os.path.dirname(__file__)
def loadpath(name):
    path = os.path.join(dirname,'save',name)
    return t.load(path,map_location=t.device('cpu'))

PARAM_ARGS = {'nin':9,'nout':1,'nh':[64]*5,'pdrop':[0.1]+[0.1]*5}


class Net_param(nn.Module):
    def __init__(
        self,
        nin,nout=1,
        nh=[128]*10,
        pdrop=[0.2]+[0.3]*10,
    ):
        super().__init__()
        fc_list = [nin, *nh, nout]
        pdrop = pdrop+[0]
        nlayer = len(nh)+1
        fc = [nn.Linear(fc_list[i],fc_list[i+1]) for i in range(nlayer)]
        act = [nn.ReLU()]*(nlayer-1) + [nn.Identity()]
        drop = [nn.Dropout(p) for p in pdrop[1:]]
        seq_list = [nn.Dropout(pdrop[0])]
        for f,a,d in zip(fc,act,drop):
            seq_list += [f,a,d]
        self.full = nn.Sequential(*seq_list)

    @t.jit.script_if_tracing
    def decode(self,output):
        return output

    @t.jit.script_if_tracing
    def forward(self,info):
        # info.shape: batch, feature
        output = self.full(info)
        return self.decode(output)


class DL_gr(Net_param):
    def __init__(self,name='dl_param-gr.pth'):
        super().__init__(**PARAM_ARGS)
        self.load_state_dict(loadpath(name))
        self.eval()

    @t.jit.script_if_tracing
    def decode(self,output):
        return 10**(output-3)

    @t.jit.script_if_tracing
    def forward(self,info):
        with t.no_grad():
            output = super().forward(info)
        return output.numpy()


class DL_h(DL_gr):
    def __init__(self,name='dl_param-h.pth'):
        super().__init__(name)

    @t.jit.script_if_tracing
    def decode(self, output):
        return output*2+1


# ============= s2s =============
S2S_ARGS = {
    'enc_args':{
        'num_inputs':10,  # 1(gamma) + 9(info)
        'num_channels':[10]*3 + [5]*3,
        'kernel_size':6
    },
    'dec_args':{
        'num_inputs':5,
        'num_channels':[10]*3 + [5]*3,
        'kernel_size':6
    }
}

class Net_s2s(nn.Module):
    def __init__(self,enc_args,dec_args):
        super().__init__()
        nmodel = 2
        self.encoder = TCN.TemporalConvNet(**enc_args)
        self.connect = nn.Linear(enc_args['num_channels'][-1],dec_args['num_inputs']-nmodel)
        self.decoder = TCN.TemporalConvNet(**dec_args)
        self.outnet = nn.Linear(dec_args['num_channels'][-1],1)

    @t.jit.script_if_tracing
    def forward(self,info,gamma,taumodel,dev='cpu'):
        # info.shape: feature
        # gamma.shape: seq
        # taumodel.shape: seq,model
        seq = len(gamma)
        gamma = gamma.reshape(1,seq,1)
        taumodel = taumodel.reshape(1,seq,2)
        info = info.reshape(1,1,-1) * t.ones(1,seq,1).to(dev)
        # print(info.shape,gamma.shape)
        x = t.cat((info,gamma),2).transpose(1,2)
        x = self.encoder(x).transpose(1,2)
        x = self.connect(x)
        x = t.cat((x,taumodel),2).transpose(1,2)
        x = self.decoder(x).transpose(1,2)
        x = self.outnet(x)
        return x


class DL_s2s(Net_s2s):
    def __init__(self,name='state_dict.pth'):
        super().__init__(**S2S_ARGS)
        self.load_state_dict(loadpath(name))
        self.eval()

    @t.jit.script_if_tracing
    def forward(self,info,gamma,taumodel):
        with t.no_grad():
            return super().forward(info,gamma,taumodel,dev='cpu')
