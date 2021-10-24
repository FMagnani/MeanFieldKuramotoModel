#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 20:38:16 2021

@author: FMagnani
GitHub: https://github.com/FMagnani
"""

import numpy as np
import matplotlib.pyplot as plt

#%%

class Integrator():
    
    def __init__(self, init_distribution, frequencies, coupling, temperature):
        

        if not (len(init_distribution) == len(frequencies)):
            raise ValueError("init_distribution and frequencies should have same length")
        
        self.phi = np.array(init_distribution)
        self.freqs = np.array(frequencies)
        self.T = temperature
        self.k = coupling
        self.coherence, self.psi = self.get_orderPars()
        
        self.scheme_dict = {
                            'euler': self.Euler,
                            'heun' : self.Heun
                            }
        
        
    def get_orderPars(self):
        """
        psi: mean phase
            psi = 1/N Sum_j phi_j
        
        r: coherence
            Amplitude (length) of the center-of-mass vector
            r = SQRT[ (1/N Sum_j cos(phi_j))^2 + (1/N Sum_j sin(phi_j))^2 ]
        """
        psi = np.mean(self.phi)
        
        cosines = np.cos(self.phi)
        sines = np.sin(self.phi)
        mean_cos = np.mean(cosines)
        mean_sin = np.mean(sines)
        coherence = np.sqrt( mean_cos**2 + mean_sin**2 )

        return coherence, psi


    def integrate(self, Dt, iterations, numerical_scheme, seed):
        """
        

        Parameters
        ----------
        Dt : float
            Time interval for integration.
        iterations : int
            Total number of iterations.
        numerical_scheme : str
            One of: 
                'euler'
                'heun'
        seed : int
            Seed for generation of random numbers for the Wiener process.

        Returns
        -------
        phi_history : np.array with shape (len(phi), iterations)
            Phase of each oscillator at each iteration.
        coherence_history : np.array with shape (len(phi))
            The coherence parameter at each iteration.
        psi_history : np.array with shape (len(phi))
            The mean phase at each iteration.
        """
                
        scheme = self.scheme_dict[numerical_scheme]
        
        coherence_history = np.array([self.coherence])
        psi_history = np.array([self.psi])
    
        # The square brackets inside are needed
        phi_history = np.array([self.phi])
         
        # For each oscillator, the whole Wiener process must be created
        # In total, we have N random walks each with length = "iterations" 
        wiener = np.random.normal(0, np.sqrt(Dt), (N,iterations))
        
        for i in range(iterations):
            
            noise = wiener[:, i]
            
            self.phi = scheme(Dt, noise)
            
            # Pay attention to the following AWFUL syntax. It's not obvious.
            phi_history = np.append(phi_history, [self.phi], axis=0)
            
            r, psi = self.get_orderPars()
            coherence_history = np.append(coherence_history, r)
            psi_history = np.append(psi_history, psi)
      
         # TODO
         # Understand the effects of 2pi modulation. I think it's different
         # for example think to a very fast oscillator and a very slow.
         # Without modulation, these two particles can get further and further
         # away, one from the other. With modulation, they can be at most 
         # far by pi. Do this affects the sinusoidal forcing between them?
         # Or maybe not since the sinusoid is inherently periodic in 2pi?
         #
         # This affects also the "speed" method that computes the distances. 
         # Correct it if needed.
         
         # A quick test showed that it doesn't affect r actually
         
        phi_history = np.mod(phi_history, 2*np.pi)
        coherence_history = np.mod(coherence_history, 2*np.pi)
        psi_history = np.mod(psi_history, 2*np.pi)
        
        return phi_history, coherence_history, psi_history

    
    # NUMERICAL SCHEMES
    
    def speed(self, phi):
        
        # Compute all combinations between phi's entries
        xx, yy = np.meshgrid(phi,phi)        
    
        # Apply all pairwise differences in parallel
        distances = xx - yy
    
        # Apply sine to all the differences in parallel
        sines = np.sin(distances)
                
        # Forcing term
        forcings = self.k*np.mean(sines, axis=1)
    
        return self.freqs + forcings

    
    def Euler(self, Dt, noise):
    
        speed = self.speed(self.phi)
    
        new_phi = self.phi + Dt*speed + self.T*noise
        
        return new_phi


    def Heun(self, Dt, noise):
        
        middle_pnt = self.Euler(Dt, noise)
        
        # Deterministic part
        deterministic_part = (self.speed(self.phi) + self.speed(middle_pnt))*.5

        new_phi = self.phi + Dt*deterministic_part + self.T*noise
        
        return new_phi
        
    ## TODO
    ## Taylor method pls
    ## I think that's a very noice integrator, this means less iterations.
    ## Less iterations means less time. Or equal iteration but with no
    ## fear to be wrong and have to re-do all from the start.
    ##
    

#%%

### SCRIPT EXAMPLE

## Number of oscillators
N = 2

## Random seed (for the Wiener processes) 
np.random.seed(123)

## Initial phases distribution
# init_phi = np.arange(0,2*np.pi,1/N) # Uniform
init_phi = np.random.normal(np.pi, 1, N) # Gaussian around pi

## Natural frequencies distribution
# freqs = np.random.normal(0,1,N) # Gaussian with 0 mean (co-rotating frame)
freqs = np.array([-5,5]) # Specific for N=2

## Coupling and Temperature parameters
K = 10
T = 0.8

system = Integrator(init_phi, freqs, K, T)

iterations = 10000

## Instead of Dt itself the resolution with wich to integrate one period is
## chosen (the smaller period is taken as a reference).
## !! BEWARE !! 
## psi should be constant during the integration. The quality of the 
## integration can be assessed in a first approximation looking at this 
time_resolution = 1000 # Number of steps used to integrate over a period
higher_freq = np.abs(freqs.max())
Dt = 1/(time_resolution*higher_freq)

phi, r, psi = system.integrate(Dt, iterations, 'heun', 7264)

## TODO
## x axis must show the length in pi units that's fundamental
## I have to understand if the behaviour is numerical (artifact) or real 
## To understand this, I could trust the model that preserves psi. That one
## tells the truth.
## If the pattern of r falling to 0 and recovering is true, then I can compute
## the average length before this to happen. If I'm good enough in math...
## And that could be an observable quantity and a prediction
## That would make me a nice guy I do believe

# fig, ax = plt.subplots()
# ax.
# plt.plot(phi, 'r-')
# ax.plot(phi[:,1], 'b-', alpha=.8)
# ax.plot(psi, 'k.')
# ax.plot(r, 'k')
plt.show()

#%%

## BIFURCATION COMPUTATION

# Instead of Dt itself the resolution with wich to integrate one period is chosen
# The smaller period is taken as a reference
time_resolution = 100 # Number of steps used to integrate over a period
higher_freq = np.abs(freqs.max())
Dt = 1/(time_resolution*higher_freq)

for (T, color) in zip([0,0.1,0.2,0.5],['k','b','g','r']):
    
    rk_relation = np.array([[0,0]])
    
    for k in np.arange(0,12, 0.4):
    
        system = Integrator(init_phi, freqs, k, T)

        iterations = 400

        _, r, _ = system.integrate(Dt, iterations, 'heun', 7264)
    
        rk_relation = np.append(rk_relation, [[k,np.mean(r[-50:])]], axis=0)

    rk_relation = rk_relation[1:,:]

    plt.plot(rk_relation[:,0], rk_relation[:,1], '-', color=color, marker='*')

plt.show()

## BIFURCATION PLOT

color='k'

plt.plot(rk_relation[:,0], rk_relation[:,1], '-', color=color, marker='*')

plt.show()


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


#%%




