#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 17:21:37 2021

@author: FMagnani
GitHub: https://github.com/FMagnani
"""

from Integrator import Integrator
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import numpy as np

#%%

"""
FIGURE none

We verify formula 45 of the report.
So we use the Cauchy distribution for the natural frequencies and we simulate
for many k. We then plot the stationary r obtained (it should be stationary...)
and such a plot should agree with the formula 45.

FAILED

"""

N = 1000 # Number of oscillators

#init_phi = np.linspace(0,2*np.pi,N) # Uniform
#init_phi = np.random.normal(np.pi, 1, N) # Gaussian around pi
init_phi = np.ones(N)

freqs = np.random.standard_cauchy(N) # g(omega) == cauchy with gamma=1
freqs -= freqs.mean()

def plot_freqs(g):
    g = g[(g>-25) & (g<25)] # Truncate distribution so it plots well
    plt.hist(g, bins=100)
    plt.show()

k_critical = 2/(np.pi) # From formula 43 recalling that gamma=1

def stable_r(k):
    if k < k_critical:
        return 0
    else:
        return np.sqrt(1-(k_critical/k))

T = 0



K = 2
system = Integrator(init_phi, freqs, K, T)

iterations = 100
Dt = 0.01
seed = 1234

phi, r, psi = system.integrate(Dt, iterations, 'taylor', seed)



time_axis = Dt*np.array(range(iterations+1))

fig, (ax1,ax2) = plt.subplots(2,1)

extract = np.random.choice(range(N), int(N/10))
ax1.plot(time_axis, phi[:,extract], linewidth=1)
ax2.plot(time_axis, r, 'k')

plt.show()

#%%

"""
FIGURE 1

Comparison betweeen schemes.
"""

N = 1000 # Number of oscillators

#init_phi = np.linspace(0,2*np.pi,N) # Uniform
init_phi = np.random.normal(np.pi, 1, N) # Gaussian around pi
#init_phi = np.ones(N)

#freqs = np.random.standard_cauchy(N) # g(omega) == cauchy with gamma=1
freqs = np.zeros(N)

T = 1
K = 4
system = Integrator(init_phi, freqs, K, T)

iterations = 100
Dt = 0.01
seed = 1234

_, r_taylor, _ = system.integrate(Dt, iterations, 'taylor', seed)
_, r_heun, _ = system.integrate(Dt, iterations, 'heun', seed)
_, r_euler, _ = system.integrate(Dt, iterations, 'euler', seed)


time_axis = Dt*np.array(range(iterations+1))

fig, ax = plt.subplots()

ax.plot(time_axis, r_taylor, 'k', label="Taylor method", alpha=.8)
ax.plot(time_axis, r_heun, 'b', label="Heun method", alpha=.8)
ax.plot(time_axis, r_euler, 'r', label="Euler method", alpha=.8)
ax.legend()

plt.show()







































