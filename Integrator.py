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
        self.w = np.array(frequencies)
        self.T = temperature
        self.k = coupling
        self.coherence, self.psi = self.get_orderPars()
        
        self.scheme_dict = {
                            'euler': self.Euler
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
                'taylor' (Taylor 3/2 is used)
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
         
        for i in range(iterations):
            
            self.phi = scheme(Dt, seed)
            
            # Pay attention to the following AWFUL syntax. It's not obvious.
            phi_history = np.append(phi_history, [self.phi], axis=0)
            
            r, psi = self.get_orderPars()
            coherence_history = np.append(coherence_history, r)
            psi_history = np.append(psi_history, psi)
            
        return phi_history, coherence_history, psi_history

    
    # NUMERICAL SCHEMES
    
    def Euler(self, Dt, seed):
    
        # Compute all combinations between phi's entries
        xx, yy = np.meshgrid(self.phi,self.phi)        
    
        # Apply all pairwise differences in parallel
        distances = xx - yy
    
        # Apply sine to all the differences in parallel
        sines = np.sin(distances)
    
        # Forcing term
        forcings = self.k*np.mean(sines, axis=1)
    
        # Noise
        np.random.seed(seed)
        noise = self.T*np.random.normal(0,np.sqrt(Dt),len(self.phi))
    
        new_phi = self.phi + Dt*(freqs + forcings + noise)
    
        # New phi modulo 2pi
        new_phi = np.mod(new_phi, 2*np.pi)
    
        return new_phi




#%%

# Script example

N = 2

np.random.seed(123)

freqs = np.random.normal(0,1,N)


init_phi_distr = np.arange(0,1,1/N) # Uniform
init_phi = np.multiply(init_phi_distr, 2*np.pi)

freqs = np.random.normal(0,1,N)

for i in range(21)[18:]:

    K = 2.2
    T = i
    system = Integrator(init_phi, freqs, K, T)
    phi, r, psi = system.integrate(0.01, 1000, 'euler', 123)
    plt.plot(phi)
    # plt.plot(r, 'k')


# system = Integrator(2, [0,1], [1,2], 3)
# phi, r, psi = system.integrate(0.01, 1000, EulerKuramoto)
#plt.plot(phi)
#plt.plot(r, 'k')





