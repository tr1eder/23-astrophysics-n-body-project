import matplotlib.pyplot as plt
import numpy as np


G = 1
M = 1

xs = np.linspace(0, 10, 100)
ys = -G*M/xs

plt.plot(xs, ys)
plt.show()
