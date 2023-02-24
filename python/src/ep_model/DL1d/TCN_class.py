import torch
import torch.nn as nn
from torch.nn.utils import weight_norm,remove_weight_norm


class Chomp1d(nn.Module):
    def __init__(self, chomp_size):
        super(Chomp1d, self).__init__()
        self.chomp_size = chomp_size

    def forward(self, x):
        return x[:, :, :-self.chomp_size].contiguous()


class TemporalBlock(nn.Module):
    def __init__(self, n_inputs, n_outputs, kernel_size, stride, dilation, padding, dropout=0.2):
        super(TemporalBlock, self).__init__()
        self.conv1 = weight_norm(nn.Conv1d(n_inputs, n_outputs, kernel_size, stride=stride, padding=padding, dilation=dilation))
        self.chomp1 = Chomp1d(padding)
        self.relu1 = nn.ReLU()
        self.dropout1 = nn.Dropout(dropout)

        self.conv2 = weight_norm(nn.Conv1d(n_outputs, n_outputs, kernel_size, stride=stride, padding=padding, dilation=dilation))
        self.chomp2 = Chomp1d(padding)
        self.relu2 = nn.ReLU()
        self.dropout2 = nn.Dropout(dropout)

        self.net = nn.Sequential(self.conv1, self.chomp1, self.relu1, self.dropout1, self.conv2, self.chomp2, self.relu2, self.dropout2)
        self.downsample = nn.Conv1d(n_inputs, n_outputs, 1) if n_inputs != n_outputs else None
        self.relu = nn.ReLU()
        self.init_weights()

    def init_weights(self):
        self.conv1.weight.data.normal_(0, 0.01)
        self.conv2.weight.data.normal_(0, 0.01)
        if self.downsample is not None:
            self.downsample.weight.data.normal_(0, 0.01)

    def remove_weight_norm(self):
        remove_weight_norm(self.conv1)
        remove_weight_norm(self.conv2)

    def forward(self, x):
        out = self.net(x)
        res = x if self.downsample is None else self.downsample(x)
        return self.relu(out + res)


class Cashed_TCN(nn.Module):
    '''`input size: (batch, input channel)`'''
    def __init__(self, num_inputs, num_channels, kernel_size=2, dropout=0):
        super(Cashed_TCN, self).__init__()
        layers = []
        inputs = []
        num_levels = len(num_channels)
        for i in range(num_levels):
            dilation_size = 2 ** i
            in_channels = num_inputs if i == 0 else num_channels[i-1]
            out_channels = num_channels[i]
            layers += [TemporalBlock(in_channels, out_channels, kernel_size, stride=1, dilation=1,padding=kernel_size-1,dropout=0.0)]
            inputs += [torch.zeros((1,in_channels,(kernel_size-1)*dilation_size+1))]

        self.network = nn.ModuleList(layers)
        self.inputs = inputs
        self.num_layer = num_levels
        self.kernel_size = kernel_size

    def remove_weight_norm(self):
        for i in range(self.num_layer):
            self.network[i].remove_weight_norm()

    def forward(self, x):
        '''
        `input size: (batch, input channel)`,
        `output size: (batch, output channel)`
        '''
        for i,net in enumerate(self.network):
            dilation_size = 2 ** i
            idx = [dilation_size*j for j in range(self.kernel_size)]
            self.inputs[i][...,:-1] = self.inputs[i][...,1:]
            self.inputs[i][...,-1] = x
            x = net(self.inputs[i][:,:,idx])[:,:,-1]
        return x