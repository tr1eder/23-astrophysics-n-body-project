import numpy as np
import matplotlib.pyplot as plt

def plotFromFile(filename):
    with open(filename) as f:
        lines = f.readlines()

    boundaries =    np.array(lines[0].split(), dtype=float)[:]
    midpoints  =    np.array(lines[1].split(), dtype=float)[:]
    computed   =    np.array(lines[2].split(), dtype=float)[:]
    rhos       =    np.array(lines[3].split(), dtype=float)[:]
    errs       =    np.array(lines[4].split(), dtype=float)[:]

    errorbars  =    lambda y : np.sqrt(y) ## one standard deviation

    plt.errorbar(midpoints, rhos, yerr=errorbars(rhos), fmt='o', label='rho w.errorbars') 
    plt.plot(midpoints, computed, label='shells')
    for x,y,err in zip(midpoints, computed, errs):
        plt.annotate('{:.1f}σ'.format(err), (x,y), textcoords="offset points", xytext=(3,5), ha='center')


    plt.yscale('log')
    # plt.ylim(0, 50)
    plt.legend()
    plt.show()


plotFromFile("plot_task1_1.txt")
