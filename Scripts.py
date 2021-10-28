#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 17:33:54 2021

@author: FMagnani
GitHub: https://github.com/FMagnani
"""


from Integrator import Integrator
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import numpy as np

#%%

### SCRIPT EXAMPLE

## Number of oscillators
N = 10

## Random seed (for the Wiener processes) 
np.random.seed(123)

## Initial phases distribution
#init_phi = np.linspace(0,2*np.pi,N) # Uniform
init_phi = np.random.normal(np.pi, 1, N) # Gaussian around pi

## Natural frequencies distribution
freqs = np.random.normal(0,1,N) # Gaussian with 0 mean (co-rotating frame)
#freqs = np.array([-5,5]) # Specific for N=2

# Actual co-rotating reference frame, since for small N real mean is not 0
freqs -= freqs.mean()

## Coupling and Temperature parameters
K = 2.5
T = 0.5

system = Integrator(init_phi, freqs, K, T)

iterations = 5000

## Instead of Dt itself the resolution with wich to integrate one period is
## chosen (the smaller period is taken as a reference).
## !! BEWARE !! 
## psi should be constant during the integration. The quality of the 
## integration can be assessed in a first approximation looking at this 
time_resolution = 100 # Number of steps used to integrate over a period
higher_freq = np.abs(freqs.max())
Dt = 1/(time_resolution*higher_freq)

phi, r, psi = system.integrate(Dt, iterations, 'heun', 7264)


modulation = 0
## Modulation of angular variables
if modulation:
    phi = np.mod(phi, 2*np.pi)  
    r = np.mod(r, 2*np.pi)
    psi = np.mod(psi, 2*np.pi)
    phi /= np.pi
    psi /= np.pi
    
#%%

## PLOTTING

time_axis = Dt*range(iterations+1)

fig, (ax1,ax2) = plt.subplots(2,1)
if modulation:
    ax1.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
    ax1.yaxis.set_major_locator(tck.MultipleLocator(base=.5))
# ax1.scatter(time_axis,phi[:,0], marker='o', c='violet', s=2)
# ax1.scatter(time_axis,phi[:,1], marker='o', c='indigo', s=2)
ax1.plot(time_axis, phi, linewidth=1)
ax2.plot(time_axis, r, 'k')

plt.show()

#%%

## Taylor method vs Heun method

## The exact trajectories are different. The amplitude of the fluctuations
## are the same. I think that these methods are to be compared in the mean
## for many noise realizations.
## NB the same noise realization has been used in the two methods, but the 
## schemes are different. Moreover, Taylor method employs two correlated
## random variables, not only a random walk. 

temperatures = [0,0.2,0.5,1]

iterations = 1000

time_resolution = 100 # Number of steps used to integrate over a period
higher_freq = np.abs(freqs.max())
Dt = 1/(time_resolution*higher_freq)

time_axis = Dt*range(iterations+1)

for T in temperatures:
   
    system = Integrator(init_phi, freqs, K, T)
    phi_h, r_h, _ = system.integrate(Dt, iterations, 'heun', 7264)
    
    system = Integrator(init_phi, freqs, K, T)
    phi_t, r_t, _ = system.integrate(Dt, iterations, 'taylor', 7264)
    
    fig, (ax1,ax2) = plt.subplots(2,1)
    ax1.plot(time_axis, phi_h, linewidth=1)
    ax2.plot(time_axis, phi_t, linewidth=1)

    fig.show()
    
    fig, (ax1,ax2) = plt.subplots(2,1)
    ax1.plot(time_axis, r_h, 'k')
    ax2.plot(time_axis, r_t, 'k')
    
    fig.show()

#%%

# POLAR PLOT

fig = plt.figure()
ax = fig.add_subplot(projection='polar')

ax.plot(psi, r, 'k')
ax.plot(psi[0], r[0], 'ro')   # Dot = Start   * --> o
ax.plot(psi[-1], r[-1], 'r*') # Star = End    * --> o

plt.show()

#%%

# ANIMATION - test

from matplotlib import animation

fig = plt.figure()
ax = plt.axes()
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])

im=plt.imshow(phi[0,:].reshape(10,10),interpolation='none')

def animate(i):
    
    im.set_array(phi[i,:].reshape(10,10))
    
    return [im]

anim = animation.FuncAnimation(fig, animate,
                               frames=iterations, interval=10, blit=True)

anim.save("name.gif", fps=30)

