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
from scipy.interpolate import interp1d

# some physical constants
cspeed = 299792.458       # nm/ps
eps_vac = 8.8541878128e-9  # pF/nm
mu_vac = 1.25663706212e-3  # ps^2/pF/nm
eV2THz = 241.79893        # 1eV = 241.79893 THz

def angle_between1(vec1, vec2):
    return np.arctan2(vec1[0]*vec2[1]-vec1[1]*vec2[0], vec1[0]*vec2[0]+vec1[1]*vec2[1])

def angle_between2(vec1, vec2):
    return np.arctan2(vec1.real*vec2.imag-vec1.imag*vec2.real, vec1.real*vec2.real+vec1.imag*vec2.imag)

def load_epsilon(path):
    data = np.loadtxt(path)
    freq = data[:,0]
    eps = interp1d(freq, data[:,1] + 1j*data[:,2])
    return freq, eps

def fill(vec, freqs):
    for i in range(0, len(vec)):
        if callable(vec[i]):
            vec[i] = vec[i](freqs)
        else:
            vec[i] = vec[i]*np.ones(len(freqs), dtype=np.complex)
    return np.array(vec)

class Material(object):
    # Material class
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

    def R1(self, theta):
        # rotation matrix for the fields
        # the fileds are defined as column vectors
        # F = (E_x, H_y, E_y, H_x)
        # The counterclockwise rotated fields are
        # F' = R1(theta) F
        #
        #             | cos  0   -sin  0   |
        # R1(theta) = | 0    cos  0    sin |
        #             | sin  0    cos  0   |
        #             | 0   -sin  0    cos |
        R = np.zeros((4,4))
        for i in range(0,4):
            R[i,i] = np.cos(theta)
        R[0,2] = -np.sin(theta)
        R[1,3] = +np.sin(theta)
        R[2,0] = +np.sin(theta)
        R[3,1] = -np.sin(theta)
        return R

    def R2(self, theta):
        # rotation matrix for the amplitudes
        # the amplitudes are defined as column vectors
        # A = (E+_x, E-_x, E+_y, E-_y)
        # The counterclockwise rotated amplitudes are
        # A' = R2(theta) A
        #
        #             | cos  0   -sin  0   |
        # R2(theta) = | 0    cos  0   -sin |
        #             | sin  0    cos  0   |
        #             | 0    sin  0    cos |
        R = np.zeros((4,4))
        for i in range(0,4):
            R[i,i] = np.cos(theta)
        R[0,2] = -np.sin(theta)
        R[1,3] = -np.sin(theta)
        R[2,0] = +np.sin(theta)
        R[3,1] = +np.sin(theta)
        return R

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
        T = np.matmul(self.a_matrix(i), self.R2(self.angle[i]))
        for k in range(i+1, f):
            T = np.matmul(self.R1(self.angle[k]-self.angle[k-1]), T)
            T = np.matmul(self.M_matrix(k), T)
        T = np.matmul(self.R1(self.angle[f]-self.angle[f-1]), T)
        T = np.matmul(self.ainv_matrix(f), T)
        T = np.matmul(self.R2(-self.angle[f]), T)
        return T

    def coeff(self, J):
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
        B[:,0] = -self.T[:,0,0]*J[0] - self.T[:,0,2]*J[1]
        B[:,1] = -self.T[:,1,0]*J[0] - self.T[:,1,2]*J[1]
        B[:,2] = -self.T[:,2,0]*J[0] - self.T[:,2,2]*J[1]
        B[:,3] = -self.T[:,3,0]*J[0] - self.T[:,3,2]*J[1]
        X = np.linalg.solve(A, B)
        AR = X[:,0:2]
        AT = X[:,2:4]
        II = np.vdot(J,J)
        IR = np.zeros((self.N), dtype = complex)
        IT = np.zeros((self.N), dtype = complex)
        for i in range(self.N):
            IR[i] = np.vdot(AR[i,:], AR[i,:])
            IT[i] = np.vdot(AT[i,:], AT[i,:])
        return AT, AR, (IT/II).real, (IR/II).real, 1 - (IT/II).real - (IR/II).real

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
        
