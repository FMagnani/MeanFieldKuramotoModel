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
    
    def __init__(self, N_oscillators, init_distribution, frequencies, coupling):
        
        # TODO
        # 1 Maybe put warning if N very large
        # 2 Phi stuff: it can be a vector or a distribution from which withdraw
        # 3 Scheme stuff: it can be a function or a name
        # 4 Noise?
        
        self.N = N_oscillators
        self.phi = np.array(init_distribution)
        self.w = np.array(frequencies)
        self.k = coupling
        self.coherence, self.psi = self.get_orderPars()
        
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

    def integrate(self, Dt, iterations, numerical_scheme):
        
        
        coherence_history = np.array([self.coherence])
        psi_history = np.array([self.psi])
    
        # The square brackets inside are needed
        phi_history = np.array([self.phi])
         
        for i in range(iterations):
            
            self.phi = numerical_scheme(self.phi, self.w, self.k, Dt)
            
            # Pay attention to the following AWFUL syntax. It's not obvious.
            phi_history = np.append(phi_history, [self.phi], axis=0)
            
            r, psi = self.get_orderPars()
            coherence_history = np.append(coherence_history, r)
            psi_history = np.append(psi_history, psi)
            
        return phi_history, coherence_history, psi_history


#%%



def EulerKuramoto(phi, freqs, coupling, Dt):

    # Compute all combinations between phi's entries
    xx, yy = np.meshgrid(phi,phi)        
    
    # Apply all pairwise differences in parallel
    distances = xx - yy
    
    # Apply sine to all the differences in parallel
    sines = np.sin(distances)
    
    # Forcing term
    forcings = coupling*np.mean(sines, axis=1)
    
    new_phi = phi + Dt*(freqs + forcings)
    
    # New phi modulo 2pi
    new_phi = np.mod(new_phi, 2*np.pi)
    
    return new_phi


#%%

# Script example

N = 100

np.random.seed(123)

freqs = np.random.normal(0,1,N)


init_phi_distr = np.arange(0,1,1/N) # Uniform
init_phi = np.multiply(init_phi_distr, 2*np.pi)

freqs = np.random.normal(0,1,N)

K = 2.2
system = Integrator(N, init_phi, freqs, K)
phi, r, psi = system.integrate(0.01, 1000, EulerKuramoto)
# plt.plot(phi)
# plt.plot(r, 'k')
# TODO
# FIX: Why r goes above 1?
# I didn't set the seed. Now r is ok. I will never fix that :(



# system = Integrator(2, [0,1], [1,2], 3)
# phi, r, psi = system.integrate(0.01, 1000, EulerKuramoto)
#plt.plot(phi)
#plt.plot(r, 'k')





