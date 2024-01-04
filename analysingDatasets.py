import matplotlib.pyplot as plt
import numpy as np

class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    def getDist(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

class Particle:
    def __init__(self, m, p, v):
        self.m = m
        self.p = p
        self.v = v

def readData(filename):
    with open(filename) as f:
        lines = f.readlines()

    particles = []
    for line in lines[1:]:
        id, m, px, py, pz, vx, vy, vz, eps = line.split()
        particles.append(Particle(float(m), Vector(float(px), float(py), float(pz)), Vector(float(vx), float(vy), float(vz))))

    return particles




if __name__ == '__main__':
    file = 'data0.txt'
    particles = readData(file)

    xs = [p.p.x for p in particles]
    ys = [p.p.y for p in particles]
    size = [p.m**(1/2)*100 for p in particles]

    plt.scatter(xs, ys, s=size, c=size, label='Particles', )
    plt.legend()
    plt.show()

