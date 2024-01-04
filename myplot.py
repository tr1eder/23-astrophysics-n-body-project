import numpy as np
import matplotlib.pyplot as plt

def plotFromFile(filename):
    with open(filename) as f:
        lines = f.readlines()

    xs =        np.array(lines[0].split(), dtype=float)[:]
    shells =    np.array(lines[1].split(), dtype=float)[:]
    rhos =      np.array(lines[2].split(), dtype=float)[:]
    errs =      np.array(lines[3].split(), dtype=float)[:]
    errorbars = lambda y : np.sqrt(y)

    plt.plot(xs, shells, label='shells')
    plt.errorbar(xs, rhos, yerr=errorbars(rhos), fmt='o', label='rho w.errorbars')
    for x,y,err in zip(xs, shells, errs):
        plt.annotate('{:.1f}σ'.format(err), (x,y), textcoords="offset points", xytext=(3,5), ha='center')


    # plt.yscale('log')
    # plt.ylim(0, 50)
    plt.legend()
    plt.show()


plotFromFile("myplot.txt")
