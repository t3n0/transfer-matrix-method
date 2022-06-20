#
# Copyright (c) 2022 Stefano Dal Forno.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from modules.utils import *

class Material(object):

    def __init__(self, N=100):
        self.N = N
        self.freqs = []
        self.epsO = []
        self.epsE = []
        self.mu = []
        self.thick = []
        self.inter = [0.0]
        self.angle = []

    def add_epsO(self, *args):
        # ordinary dielectric function
        for arg in args:
            if isinstance(arg, str):
                f, e = load_epsilon(arg)
                self.freqs.append(f)
                self.epsO.append(e)
            else:
                self.epsO.append(arg)
        
    def add_epsE(self, *args):
        # extraordinary dielectric function
        for arg in args:
            if isinstance(arg, str):
                f, e = load_epsilon(arg)
                self.freqs.append(f)
                self.epsE.append(e)
            else:
                self.epsE.append(arg)

    def add_mu(self, *args):
        for arg in args:
            self.mu.append(arg)

    def add_layer(self, *args):
        for arg in args:
            self.thick.append(arg)

    def add_angle(self, *args):
        for arg in args:
            self.angle.append(arg)

    def M_matrix(self, i):
        M = np.zeros((self.N, 4, 4), dtype=np.complex)
        kzO   = self.freqs*np.sqrt(self.epsO[i]*self.mu[i])*self.thick[i-1]
        kzE   = self.freqs*np.sqrt(self.epsE[i]*self.mu[i])*self.thick[i-1]
        betaO = np.sqrt(self.epsO[i]/self.mu[i])
        betaE = np.sqrt(self.epsE[i]/self.mu[i])
        M[:, 0, 0] = np.cos(kzO)
        M[:, 0, 1] = 1j*np.sin(kzO)/betaO
        M[:, 1, 0] = 1j*np.sin(kzO)*betaO
        M[:, 1, 1] = np.cos(kzO)
        M[:, 2, 2] = np.cos(kzE)
        M[:, 2, 3] = -1j*np.sin(kzE)/betaE
        M[:, 3, 2] = -1j*np.sin(kzE)*betaE
        M[:, 3, 3] = np.cos(kzE)
        return M

    def a_matrix(self, i):
        a = np.zeros((self.N, 4, 4), dtype=np.complex)
        kzO   = 1j*self.freqs*np.sqrt(self.epsO[i]*self.mu[i])*self.inter[i]
        kzE   = 1j*self.freqs*np.sqrt(self.epsE[i]*self.mu[i])*self.inter[i]
        betaO = np.sqrt(self.epsO[i]/self.mu[i])
        betaE = np.sqrt(self.epsE[i]/self.mu[i])
        a[:, 0, 0] = np.exp(kzO)
        a[:, 0, 1] = np.exp(-kzO)
        a[:, 1, 0] = np.exp(kzO)*betaO
        a[:, 1, 1] = -np.exp(-kzO)*betaO
        a[:, 2, 2] = np.exp(kzE)
        a[:, 2, 3] = np.exp(-kzE)
        a[:, 3, 2] = -np.exp(kzE)*betaE
        a[:, 3, 3] = np.exp(-kzE)*betaE
        return a

    def ainv_matrix(self, i):
        ainv = np.zeros((self.N, 4, 4), dtype=np.complex)
        kzO   = 1j*self.freqs*np.sqrt(self.epsO[i]*self.mu[i])*self.inter[i-1]
        kzE   = 1j*self.freqs*np.sqrt(self.epsE[i]*self.mu[i])*self.inter[i-1]
        betaO = np.sqrt(self.epsO[i]/self.mu[i])
        betaE = np.sqrt(self.epsE[i]/self.mu[i])
        ainv[:, 0, 0] = np.exp(-kzO)/2
        ainv[:, 0, 1] = np.exp(-kzO)/2/betaO
        ainv[:, 1, 0] = np.exp(kzO)/2
        ainv[:, 1, 1] = -np.exp(kzO)/2/betaO
        ainv[:, 2, 2] = np.exp(-kzE)/2
        ainv[:, 2, 3] = -np.exp(-kzE)/2/betaE
        ainv[:, 3, 2] = np.exp(kzE)/2
        ainv[:, 3, 3] = np.exp(kzE)/2/betaE
        return ainv

    def T_matrix(self, i, f):
        T = np.matmul(self.a_matrix(i), R2(self.angle[i]))
        for k in range(i+1, f):
            T = np.matmul(R1(self.angle[k]-self.angle[k-1]), T)
            T = np.matmul(self.M_matrix(k), T)
        T = np.matmul(R1(self.angle[f]-self.angle[f-1]), T)
        T = np.matmul(self.ainv_matrix(f), T)
        T = np.matmul(R2(-self.angle[f]), T)
        return T

    def coefficients(self, IA):
        # Calculates the Transmission, Reflection and Absorption coefficients
        # given the IA : complex incident amplitudes (Jones vector)
        # input
        #     IA : ndarray[Ex, EY]         constant 2D complex array
        #     IA : ndarray[self.N, Ex, EY] frequency dependent 2D complex array
        # returns
        #     complex electric field amplitudes: IA, TA and RA, incident, transmitted and reflected amplitudes
        #     (real) electric field intensities: TI, RI and AI, transmitted, reflected and absorbed intensities

        if IA.shape == (2,):
            IA = np.ones((self.N,2))*IA
        II = np.zeros((self.N), dtype = complex) # incident intensity
        RI = np.zeros((self.N), dtype = complex) # reflected intensity
        TI = np.zeros((self.N), dtype = complex) # transmitted intensity
        A = np.zeros((self.N, 4, 4), dtype=np.complex)
        B = np.zeros((self.N, 4), dtype=np.complex)
        A[:,0,0] = self.T[:,0,1]
        A[:,0,1] = self.T[:,0,3]
        A[:,0,2] = -1
        A[:,0,3] = 0
        A[:,1,0] = self.T[:,1,1]
        A[:,1,1] = self.T[:,1,3]
        A[:,1,2] = 0
        A[:,1,3] = 0
        A[:,2,0] = self.T[:,2,1]
        A[:,2,1] = self.T[:,2,3]
        A[:,2,2] = 0
        A[:,2,3] = -1
        A[:,3,0] = self.T[:,3,1]
        A[:,3,1] = self.T[:,3,3]
        A[:,3,2] = 0
        A[:,3,3] = 0
        B[:,0] = -self.T[:,0,0]*IA[:,0] - self.T[:,0,2]*IA[:,1]
        B[:,1] = -self.T[:,1,0]*IA[:,0] - self.T[:,1,2]*IA[:,1]
        B[:,2] = -self.T[:,2,0]*IA[:,0] - self.T[:,2,2]*IA[:,1]
        B[:,3] = -self.T[:,3,0]*IA[:,0] - self.T[:,3,2]*IA[:,1]
        X = np.linalg.solve(A, B)
        RA = X[:,0:2]
        TA = X[:,2:4]
        for i in range(self.N):
            II[i] = np.vdot(IA[i,:], IA[i,:])
            RI[i] = np.vdot(RA[i,:], RA[i,:])
            TI[i] = np.vdot(TA[i,:], TA[i,:])
        return IA, TA, RA, (TI/II).real, (RI/II).real, 1 - (TI/II).real - (RI/II).real

    def TMM(self, i, f):
        for j in range(len(self.thick)):
            self.inter.append(self.inter[j]+self.thick[j])
        self.freqs = np.array(self.freqs)
        freqmin = max(self.freqs[:,0])
        freqmax = min(self.freqs[:,-1])
        self.freqs = np.linspace(freqmin, freqmax, self.N)
        self.epsO = eps_vac*fill(self.epsO, self.freqs)
        self.epsE = eps_vac*fill(self.epsE, self.freqs)
        self.mu = mu_vac*fill(self.mu, self.freqs)
        self.freqs = self.freqs * eV2THz
        self.T = self.T_matrix(i, f)
        
