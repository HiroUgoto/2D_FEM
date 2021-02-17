import numpy as np
import scipy.signal
import matplotlib.pyplot as plt

fsamp = 10000
fp = 0.196
nstart = int(2.0/fp*fsamp)
ndis = int(0.5/fp*fsamp)

vs_list = ["00", "01", "02", "03", "04", "05", "06", "08", "10"]

for vs in vs_list:
    _,disp = np.loadtxt("z0_vs" + vs + "_md.disp",usecols=(0,1),unpack=True)
    peak_p,_ = scipy.signal.find_peaks(disp[nstart:], distance=ndis)
    peak_m,_ = scipy.signal.find_peaks(-disp[nstart:], distance=ndis)
    peaks = np.hstack([peak_p,peak_m])
    max_disp = np.mean(np.abs(disp[nstart+peaks]))
    print(int(vs),max_disp)

    # plt.plot(disp)
    # plt.plot(nstart+peaks,disp[nstart+peaks],"x")
    # plt.plot(np.ones_like(disp)*max_disp)
    # plt.plot(-np.ones_like(disp)*max_disp)
    # plt.show()

print("")

for vs in vs_list:
    _,disp = np.loadtxt("z0_vs" + vs + ".disp",usecols=(0,1),unpack=True)
    peak_p,_ = scipy.signal.find_peaks(disp[nstart:], distance=ndis)
    peak_m,_ = scipy.signal.find_peaks(-disp[nstart:], distance=ndis)
    peaks = np.hstack([peak_p,peak_m])
    max_disp = np.mean(np.abs(disp[nstart+peaks]))
    print(int(vs),max_disp)

    plt.plot(disp)
    plt.plot(nstart+peaks,disp[nstart+peaks],"x")
    plt.plot(np.ones_like(disp)*max_disp)
    plt.plot(-np.ones_like(disp)*max_disp)
    plt.show()
