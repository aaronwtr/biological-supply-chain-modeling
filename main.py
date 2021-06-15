import numpy as np
import matplotlib.pyplot as plt
import time
import utils.output_handler as oh
from numba import njit

# This code is exclusively for a 1D simulation


start = time.time()


@njit()
def petri_rod(Q, N, dt, dur, samples):

    # #--------------parameters---------------------- 

    DBmax = 5.25 * 10 ** 5  # um^2/h
    Kv = .001  # um^-2
    Dp = 1  # um^2/h
    eta = 8 * 10 ** 4  # um^2/h
    beta = 80  # 1
    kl = 2  # 00          # um^2/h
    gmax = 6  # 1/h
    Kn = .1  # um^-2
    Kb = .1  # um^-2
    n0 = 1  # um^-2
    Length_petri = 1000  # um
    B_0 = 1  # um^-2
    L_0 = 1  # um^-2
    P_0 = 1  # um^-2

    # --------------Condensed constants----------------------

    Drate = DBmax * (n0 / (n0 + Kv)) ** 2  # um^2/h
    Dpr = Dp / Drate  # 1
    G = gmax * n0 / (n0 + Kn) / Drate  # um^-2
    Omega = kl * n0 / Drate  # um^-2
    h = eta * Kb / Drate  # um^-2
    Theta = beta * kl * n0 * eta * Kb / Drate ** 2  # um^-4

    q = np.zeros((3, N))
    q[0, :] = np.log(Q[0, :] / B_0)
    q[1, :] = np.log(Q[1, :] / L_0)
    q[2, :] = np.log(Q[2, :] * h / P_0)

    Results = np.zeros((3, N, samples + 1))

    Results[:, :, 0] = q[:3, :]

    Steps = np.round(dur / (dt * samples))
    DX = 1 / (4 * (Length_petri / N) ** 2)
    Dtau = Drate * dt
    DtauG = Dtau * G
    DtauDX = DX * Dtau
    DtauDXDpr = DtauDX * Dpr
    Dtauh = Dtau * h
    DtauOmega = Dtau * Omega
    logDtauTheta = np.log(Dtau * Theta)
    hKb = h * Kb
    operator = np.zeros((3, N))

    for i in range(samples):
        for j in range(Steps):
            Dtau_over_ExpB_plus_ExpL_plus_Kb = Dtau / (np.exp(q[0, :]) + np.exp(q[1, :]) + Kb)

            operator[:, 0] = (q[:, 1] - q[:, 0] + 2) ** 2 - 4
            operator[:, 1:N - 1] = (q[:, 2:N] - q[:, 0:N - 2] + 2) ** 2 + 8 * (q[:, 0:N - 2] - q[:, 1:N - 1]) - 4
            operator[:, N - 1] = (q[:, N - 2] - q[:, N - 1] + 2) ** 2 - 4

            b1 = q[0, :] + DtauG - np.exp(q[2, :]) * Dtau_over_ExpB_plus_ExpL_plus_Kb + operator[0, :] * DtauDX
            l1 = q[1, :] + np.exp(
                q[2, :] + q[0, :] - q[1, :]) * Dtau_over_ExpB_plus_ExpL_plus_Kb - DtauOmega + operator[1, :] * DtauDX
            q[2, :] += -Dtauh + hKb * Dtau_over_ExpB_plus_ExpL_plus_Kb + np.exp(
                q[1, :] - q[2, :] + logDtauTheta) + operator[2, :] * DtauDXDpr

            q[0, :] = b1
            q[1, :] = l1

        Results[:, :, i + 1] = q

    return B_0 * np.exp(Results[0, :, :]), L_0 * np.exp(Results[1, :, :]), P_0 * np.exp(Results[2, :, :]) / h,


# ------------------Simulation Variables---------------------

N = 50  # discretisaion in X (N)
dt = 10 ** (-8)  # time steps in h
dur = .30  # duration of total simulation in h
samples = 300  # how often one wants data to be returned

# parameters for bacteria, infected bacteria and phages in # um^-2

B = np.ones((N)) * 10 ** (-6);
L = np.ones((N)) * 10 ** (-6);
P = np.ones((N)) * 10 ** (-6)  # **(-6) = 0
B[0] = 100
P[0] = 10 ** (-2)

# convert BLP to blp

Q = np.zeros((3, N))

Q[0, :] = B
Q[1, :] = L
Q[2, :] = P

Bresults, Lresults, Presults = petri_rod(Q, N, dt, dur, samples)

sr = True
if sr:
    oh.save_results(Bresults, Lresults, Presults)

print("The simulation took " + str(time.time() - start) + " seconds")
