import numpy as np
import matplotlib.pyplot as plt
import ast

def plotFromFile(filename):
    with open(filename) as f:
        lines = f.readlines()

    boundaries =    np.array(lines[0].split(), dtype=float)[:]
    analytical =    np.array(lines[1].split(), dtype=float)[:]
    computed   =    np.array(lines[2].split(), dtype=float)[:]
    angleOff   =    np.array(lines[3].split(), dtype=float)[:]
    scatter    =    list(map(list, zip(*ast.literal_eval(lines[4]))))

    # print (scatter)

    fig, ax1 = plt.subplots()

    ax1.plot(boundaries, analytical, label='Analytical Newton II')
    ax1.plot(boundaries, computed, label='Direct forces')
    ax1.scatter(scatter[0], scatter[1], s=3, c='k', label='Scatter points')
    # ax1.set_yscale('symlog')

    ax2 = ax1.twinx()
    ax2.plot(boundaries, np.degrees(angleOff), 'r--', label='Angle offset (rad)')
    ax2.tick_params('y', colors='r')
    ax2.set_ylim(0, 10)
    ax2.yaxis.set_major_formatter('{x}°')

    ax1.legend()
    ax2.legend()
    plt.show()


plotFromFile("plot_task1_2.txt")
