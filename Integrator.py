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

# Made in general reference frame -> d/dt (Mean phase) is not 0
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
         
        for i in range(iterations-2):
            
            self.phi = scheme(Dt, seed)
            
            # Pay attention to the following AWFUL syntax. It's not obvious.
            phi_history = np.append(phi_history, [self.phi], axis=0)
            
            r, psi = self.get_orderPars()
            coherence_history = np.append(coherence_history, r)
            psi_history = np.append(psi_history, psi)
            
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

    
    def Euler(self, Dt, seed):
    
        # Noise
        np.random.seed(seed)
        noise = self.T*np.random.normal(0,np.sqrt(Dt),len(self.phi))
    
        speed = self.speed(self.phi)
    
        new_phi = self.phi + Dt*speed + noise
        
        return new_phi


    def Heun(self, Dt, seed):
        
        middle_pnt = self.Euler(Dt, seed)
        
        # Deterministic part
        deterministic_part = (self.speed(self.phi) + self.speed(middle_pnt))*.5

        # Noise
        np.random.seed(seed)
        noise = self.T*np.random.normal(0,np.sqrt(Dt),len(self.phi))
        
        new_phi = self.phi + Dt*deterministic_part + noise
        
        return new_phi
        

#%%

# Script example

# Number of oscillators
N = 10

np.random.seed(123)

# Initial phases distribution
# init_phi_distr = np.arange(0,1,1/N) # Uniform
# init_phi = np.multiply(init_phi_distr, 2*np.pi)

init_phi = np.random.normal(np.pi, 1, N)

# Natural frequencies distribution
freqs = np.random.normal(10,1,N)

# System parameters
K = 10
T = 0

system = Integrator(init_phi, freqs, K, T)


iterations = 2000

# Instead of Dt itself the resolution with wich to integrate one period is chosen
time_resolution = 80 # Number of steps used to integrate over a period
average_freq = np.abs(freqs.mean())
Dt = 1/(time_resolution*average_freq)

phi, r, psi = system.integrate(Dt, iterations, 'heun', 123)
# plt.plot(phi, 'r+')
# plt.plot(psi, 'k+')
# plt.plot(r, 'b')

#%%

fig = plt.figure()
ax = fig.add_subplot(projection='polar')

#ax.plot(phi[:,0], speed[:,0], 'r')
#ax.plot(phi[:,1], speed[:,1], 'b')
#ax.plot(phi[:,9], speed[:,9], 'k')

ax.plot(psi, r, 'r')
