import math
import numpy as np

# G = 6.67*10**-11  # Gravitational constant (set to 1 for simplicity)
# M_solar_mass = 2*10**30  # Solar mass in arbitrary units
# AU = 1.5*10**11  # Astronomical unit in arbitrary units

# # Calculate the orbital period T using the formula
# T = 2 * math.pi * math.sqrt(AU**3 / (G * M_solar_mass))

# print(f"Orbital period T: {T} time units")
# print (T/60/60/24/365)

d = 1.5*10**11
mass = 5.9722*10**24
# mass = 2*10**30


t = np.sqrt(d**3 / (6.67 * 10**-11 * mass))
v = d/t
print (t)
print (v)
print (2.978*10**4 / v)

print (t/60/60/24/365)

# print (2*np.pi*t)

# solMass = 1.98892 * 10**30
# earthMass = 5.9722 * 10**24

# print (solMass/earthMass)



"""
    6.67*10^-11 m^3/(kg*s^2) = 1 * [d^3]/([m]*[t]^2)

    eg. 6.67*(10**-11) / (d**3) * mass = 1/t**2

    t = np.sqrt(d**3 / (6.67 * 10**-11 * mass))




"""