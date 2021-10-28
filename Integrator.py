#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 20:38:16 2021

@author: FMagnani
GitHub: https://github.com/FMagnani
"""

import numpy as np

#%%

class Integrator():
    
    def __init__(self, init_distribution, frequencies, coupling, temperature):
        

        if not (len(init_distribution) == len(frequencies)):
            raise ValueError("init_distribution and frequencies should have same length")
        
        self.N = np.array(init_distribution).shape[0]
        self.phi = np.array(init_distribution)
        self.freqs = np.array(frequencies)
        self.T = temperature
        self.k = coupling
        self.coherence, self.psi = self.get_orderPars()
        
        self.scheme_dict = {
                            'euler': self.Euler,
                            'heun': self.Heun,
                            'taylor': self.Taylor
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

        N = self.N


        ## TODO
        ## Update the following description
        
        ## The noise realizations must be created for all the iterations 
        ## before the loop. If we create a single one at each loop iteration,
        ## since the seed has been chosen, they will be all equal.
        ## Furthermore, some schemes as the Taylor scheme needs two random
        ## variables correlated. They are created via a multivariate gaussian
        ## distribution. 
        ## Later in the scheme function, Euler and Heun methods will employ
        ## only the Wiener part (the firs of the two correlated variables)
        ## while Taylor both.
        
        means = [0,0]
        covs = [[np.sqrt(Dt),     0.5*(Dt**2)],
                [0.5*(Dt**2),     (Dt**3)/3  ]]
        noise_vars = [np.random.multivariate_normal(means, covs, N) for i in range(iterations)]
            
        for i in range(iterations):

            noise = noise_vars[i]
            
            self.phi = scheme(Dt, noise)
            
            # Pay attention to the following AWFUL syntax. It's not obvious.
            phi_history = np.append(phi_history, [self.phi], axis=0)
            
            r, psi = self.get_orderPars()
            coherence_history = np.append(coherence_history, r)
            psi_history = np.append(psi_history, psi)
         
        # phi_history = np.mod(phi_history, 2*np.pi)
        # coherence_history = np.mod(coherence_history, 2*np.pi)
        # psi_history = np.mod(psi_history, 2*np.pi)
        
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
        
        # Noise comes with Wiener part + first order. 
        # Only Wiener part is used here.
        # noise is a vector with a different perturbation for each oscillator.
        wiener = noise[:,0]
        
        speed = self.speed(self.phi)
        
        new_phi = self.phi + Dt*speed + self.T*wiener
        
        return new_phi


    def Heun(self, Dt, noise):

        # Noise comes with Wiener part + first order. 
        # Only Wiener part is used here.
        # noise is a vector with a different perturbation for each oscillator.
        wiener = noise[:,0]
        
        middle_pnt = self.Euler(Dt, noise)
        
        # Deterministic part
        deterministic_part = (self.speed(self.phi) + self.speed(middle_pnt))*.5

        new_phi = self.phi + Dt*deterministic_part + self.T*wiener
        
        return new_phi
        
    
    ## TODO
    ## Assess if Taylor method is right
    
    def speed_prime(self, phi):
        
        # Compute all combinations between phi's entries
        xx, yy = np.meshgrid(phi,phi)        
    
        # Apply all pairwise differences in parallel
        distances = xx - yy
    
        # Apply sine to all the differences in parallel
        cosines = np.cos(distances)
                
        # Forcing term
        forcings = -self.k*np.mean(cosines, axis=1)
    
        return forcings    
    
    
    def Taylor(self, Dt, noise):
        
        noise = np.array(noise)
        
        DW = noise[:,0]
        DZ = noise[:,1]

        f = self.speed(self.phi)
        f_prime = self.speed_prime(self.phi)
        
        new_phi = self.phi + f*Dt + self.T*DW +\
                  0.5*Dt*Dt*f*( f_prime + 0.5*self.T*self.T ) +\
                  self.T*f_prime*DZ

        return new_phi


#%%




