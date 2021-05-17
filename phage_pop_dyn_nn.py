from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

import deepxde as dde
from deepxde.backend import tf

def main():
    def pde(x, y):
        b = 100     # initial number of bacteria
        p = 10      # initial number of phages
        l = 0       # initial number of infected bacteria

        T = b + l
        F = p / (T + Kb)

        db_tau = dde.grad.jacobian(y, x, i=0, j=1)
        dl_tau = dde.grad.jacobian(y, x, i=0, j=1)
        dp_tau = dde.grad.jacobian(y, x, i=0, j=1)
        db_xx = dde.grad.hessian(y, x, i=0, j=0)
        dl_xx = dde.grad.hessian(y, x, i=0, j=0)
        dp_xx = dde.grad.hessian(y, x, i=0, j=0)
        return [db_tau - (db_xx + (G - F)*b),
                dl_tau - (dl_xx - Omega*l + b*F),
                dp_tau - (Dpr*dp_xx - T*h*F + l*Theta)]


    def func(x):
        return 0

    geom = dde.geometry.Interval(-10, 10)
    timedomain = dde.geometry.TimeDomain(0, 30)
    geomtime = dde.geometry.GeometryXTime(geom, timedomain)

    bc = dde.DirichletBC(geomtime, func, lambda _, on_boundary: on_boundary)
    ic = dde.IC(geomtime, func, lambda _, on_initial: on_initial)
    data = dde.data.TimePDE(
        geomtime,
        pde,
        [bc, ic],
        num_domain=40,
        num_boundary=20,
        num_initial=10,
        solution=func,
        num_test=10000,
    )

    layer_size = [2] + [32] * 3 + [1]
    activation = "tanh"
    initializer = "Glorot uniform"
    net = dde.maps.FNN(layer_size, activation, initializer)

    model = dde.Model(data, net)

    model.compile("adam", lr=0.001, metrics=["l2 relative error"])
    losshistory, train_state = model.train(epochs=10000)

    dde.saveplot(losshistory, train_state, issave=True, isplot=True)


if __name__ == "__main__":
    '''
    Parameters
    '''
    DBmax = 5.25 * 10 ** 5  # um^2/h
    Kv = .001  # um^-2
    Kc = .02  # um^-2
    alphac = 2  # 1
    Dp = 1  # um^2/h
    Dn = 4.5 * 10 ** 6  # um^2/h
    eta = 8 * 10 ** 4  # um^2/h
    beta = 80  # 1
    kl = 2  # 00             # um^2/h
    gmax = 6  # 1/h
    Kb = .1  # um^-2
    Lambda = .2  # 1
    Kb = .1  # um^-2
    B0 = 10 ** 7  # 1
    P0 = 1.5 * 10 ** 8  # 1
    n0 = 1  # um^-2
    Length_petri = 1000  # um
    B_0 = 1  # um^-2
    L_0 = 1  # um^-2
    P_0 = 1  # um^-2
    B_min = 10 ** (-10)  # um^-2 = 1 cm^-2
    L_min = 10 ** (-10)  # um^-2 = 1 cm^-2
    P_min = 10 ** (-10)  # um^-2 = 1 cm^-2

    '''
    Condensed constants
    '''
    Drate = DBmax * (n0 / (n0 + Kv)) ** 2  # um^2/h
    Dpr = Dp / Drate  # 1
    G = gmax * n0 / (n0 + Kb) / Drate  # um^-2
    Omega = kl * n0 / Drate  # um^-2
    h = eta * Kb / Drate  # um^-2
    Theta = beta * kl * n0 * eta * Kb / Drate ** 2  # um^-4

    # '''
    # Vector of constants
    # '''
    # Vc = [Drate, Dpr, G, Omega, h, Theta, Length_petri]
    #
    #
    # '''
    # Simulation variables
    # '''
    # N = 10  # discretisaion in X (N)
    # dt = 10 ** (-8)  # time steps in h
    # dur = .0030  # duration of total simulation in h
    # samples = 30  # how often one wants data to be returned
    #
    # Bresults = np.zeros((N, samples + 1))
    # Lresults = np.zeros((N, samples + 1))
    # Presults = np.zeros((N, samples + 1))
    #
    #
    # '''
    # Parameters for bacteria, infected bacteria and phages in micrometer^(-2)
    # '''
    # min_val = np.finfo(float).eps
    # B = min_val * np.ones(N)
    # L = min_val * np.ones(N)
    # P = min_val * np.ones(N)
    # B = np.ones(N)
    # L = np.ones(N) * 10 ** (-6)
    # P = np.ones(N) * 10 ** (-2)
    # B[0] = 100

    main()



