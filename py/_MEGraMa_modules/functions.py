import numpy as np

def rot_mat(Ang,N):

    R = np.array([[np.cos(Ang) + N[0]*N[0]*(1.0-np.cos(Ang)), N[0]*N[1]*(1.0-np.cos(Ang)) - N[2]*np.sin(Ang), N[0]*N[2]*(1.0-np.cos(Ang))+N[1]*np.sin(Ang)],
    [N[1]*N[0]*(1.0-np.cos(Ang)) + N[2]*np.sin(Ang), np.cos(Ang) + N[1]*N[1]*(1.0-np.cos(Ang)), N[1]*N[2]*(1.0-np.cos(Ang))-N[0]*np.sin(Ang)],
    [N[2]*N[0]*(1.0-np.cos(Ang)) - N[1]*np.sin(Ang), N[2]*N[1]*(1.0-np.cos(Ang)) + N[0]*np.sin(Ang), np.cos(Ang) + N[2]*N[2]*(1.0-np.cos(Ang))]])

    return R


def exp_decay_fit(x, a, b, c):

    return a * np.exp(-b * x) + c


def log_fit(x, a, b, c):

    return a * np.log(b * x) + c
