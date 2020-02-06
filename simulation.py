import numpy as np
from scipy.integrate import ode

from model.name2idx import parameters as C
from model.name2idx import variables as V
from model.param_const import f_params
from model.initial_condition import initial_values
from model.differential_equation import diffeq


def solveode(diffeq, y0, tspan, args):
    sol = ode(diffeq)
    sol.set_integrator(
        'vode', method='bdf', min_step=1e-8, with_jacobian=True
    )
    sol.set_initial_value(y0, tspan[0])
    sol.set_f_params(args)

    T = [tspan[0]]
    Y = [y0]

    while sol.successful() and sol.t < tspan[-1]:
        sol.integrate(sol.t+1.)
        T.append(sol.t)
        Y.append(sol.y)

    return np.array(T), np.array(Y)


def get_steady_state(diffeq, y0, tspan, args,
                     steady_state_time=10000, steady_state_eps=1e-8):
    sol = ode(diffeq)
    sol.set_integrator(
        'vode', method='bdf', min_step=1e-8, with_jacobian=True
    )
    sol.set_initial_value(y0, tspan[0])
    sol.set_f_params(args)

    T = [tspan[0]]
    Y = [y0]

    while sol.successful() and sol.t < steady_state_time:
        sol.integrate(steady_state_time, step=True)
        if tspan[-1] < sol.t and np.all(np.abs(sol.y - Y[-1]) < steady_state_eps):
            break
        else:
            T.append(sol.t)
            Y.append(sol.y)

    return T[-1], Y[-1]


class Simulation(object):

    tspan = range(5401)  # Unit time: 1 sec.
    condition = 2

    t = np.array(tspan)/60.  # sec. -> min. (plot_func.py)

    PMEK_cyt = np.empty((len(tspan), condition))
    PERK_cyt = np.empty((len(tspan), condition))
    PRSK_wcl = np.empty((len(tspan), condition))
    PCREB_wcl = np.empty((len(tspan), condition))
    DUSPmRNA = np.empty((len(tspan), condition))
    cFosmRNA = np.empty((len(tspan), condition))
    cFosPro = np.empty((len(tspan), condition))
    PcFos = np.empty((len(tspan), condition))

    x = f_params()
    y0 = initial_values()

    # get steady state -- preprocess
    x[C.Ligand] = x[C.no_ligand]
    (T_steady_state, Y_steady_state) = get_steady_state(diffeq, y0, tspan, tuple(x))
    if T_steady_state < tspan[-1]:
        print('Simulation failed.')
    else:
        y0 = Y_steady_state[:]
    # add ligand
    for i in range(condition):
        if i == 0:
            x[C.Ligand] = x[C.EGF]
        elif i == 1:
            x[C.Ligand] = x[C.HRG]

        (T, Y) = solveode(diffeq, y0, tspan, tuple(x))

        if T[-1] < tspan[-1]:
            print('Simulation failed.')
        else:
            PMEK_cyt[:, i] = Y[:, V.ppMEKc]
            PERK_cyt[:, i] = Y[:, V.pERKc] + Y[:, V.ppERKc]
            PRSK_wcl[:, i] = Y[:, V.pRSKc] + Y[:, V.pRSKn]*(x[C.Vn]/x[C.Vc])
            PCREB_wcl[:, i] = Y[:, V.pCREBn]*(x[C.Vn]/x[C.Vc])
            DUSPmRNA[:, i] = Y[:, V.duspmRNAc]
            cFosmRNA[:, i] = Y[:, V.cfosmRNAc]
            cFosPro[:, i] = (Y[:, V.pcFOSn] + Y[:, V.cFOSn]) * (x[C.Vn]/x[C.Vc]) \
                            + Y[:, V.cFOSc] + Y[:, V.pcFOSc]
            PcFos[:, i] = Y[:, V.pcFOSn]*(x[C.Vn]/x[C.Vc]) + Y[:, V.pcFOSc]
