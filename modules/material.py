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
import os
from modules.utils import *


class Material(object):

    def __init__(self, materials={'air': 1.0}, efields={'x_pol': 1.0}, default=[0.0, 10.0, 1000]):
        self.freqs = freqMesh(materials, efields, default)
        self.N = len(self.freqs)
        self.materials = {}
        self.IA = {}
        self.T = {}
        self.coeff = {}

        WORK_DIR = os.getcwd()
        # Create dielectric response functions for each material
        units = np.array([eps_vac, eps_vac, mu_vac])
        for mat in materials:
            if isinstance(materials[mat], str):
                _, epsO = load_file(os.path.join(
                    WORK_DIR, materials[mat], 'epsO.txt'))
                _, epsE = load_file(os.path.join(
                    WORK_DIR, materials[mat], 'epsE.txt'))
                self.materials[mat] = [epsO, epsE, 1.0]
            elif isinstance(materials[mat], int) or isinstance(materials[mat], float):
                self.materials[mat] = [materials[mat], materials[mat], 1.0]
            elif isinstance(materials[mat], list):
                eps = materials[mat][0] + 1j*materials[mat][1]
                self.materials[mat] = [eps, eps, 1.0]
        for mat in self.materials:
            self.materials[mat] = units[:, np.newaxis] * \
                fill(self.materials[mat], self.freqs)

        # Create incoming electric field amplitudes
        x_pol = np.array([1, 0],   dtype=complex)
        y_pol = np.array([0, 1],   dtype=complex)
        left_pol = np.array([1, -1j], dtype=complex)
        right_pol = np.array([1, 1j],  dtype=complex)
        for ef in efields:
            amp = efields[ef]
            if type(amp) == int or type(amp) == float:
                self.IA[ef] = amp * x_pol * np.ones((self.N, 2), dtype=complex)
            elif isinstance(amp, list):
                if amp[1] == 'x':
                    self.IA[ef] = amp[0] * x_pol * \
                        np.ones((self.N, 2), dtype=complex)
                elif amp[1] == 'y':
                    self.IA[ef] = amp[0] * y_pol * \
                        np.ones((self.N, 2), dtype=complex)
                elif amp[1] == 'left':
                    self.IA[ef] = amp[0] * left_pol * \
                        np.ones((self.N, 2), dtype=complex)
                elif amp[1] == 'right':
                    self.IA[ef] = amp[0] * right_pol * \
                        np.ones((self.N, 2), dtype=complex)
                else:
                    raise NotImplementedError('Polarisation not defined.')
            elif isinstance(amp, str):
                _, Ex = load_file(os.path.join(WORK_DIR, amp, 'Ex.txt'))
                _, Ey = load_file(os.path.join(WORK_DIR, amp, 'Ey.txt'))
                self.IA[ef] = np.zeros((self.N, 2), dtype=complex)
                self.IA[ef][:, 0] = Ex(self.freqs)
                self.IA[ef][:, 1] = Ey(self.freqs)

    def setGeometry(self, layers, thicks, angles):
        self.epsO = []
        self.epsE = []
        self.mu = []
        self.thick = []
        self.inter = [0.0]
        self.angle = []
        assert len(layers) == (len(thicks) + 2) == len(angles)
        self.epsO.append(self.materials[layers[0]][0])
        self.epsE.append(self.materials[layers[0]][1])
        self.mu.append(self.materials[layers[0]][2])
        self.angle.append(angles[0])
        for i in range(1, len(layers)-1):
            if (isinstance(thicks[i-1], int) or isinstance(thicks[i-1], float)) and (isinstance(angles[i], int) or isinstance(angles[i], float)):
                self.epsO.append(self.materials[layers[i]][0])
                self.epsE.append(self.materials[layers[i]][1])
                self.mu.append(self.materials[layers[i]][2])
                self.angle.append(angles[i])
                self.thick.append(thicks[i-1])
                self.inter.append(self.inter[-1]+self.thick[i-1])
            elif isinstance(thicks[i-1], list) and isinstance(angles[i], list):
                thick_unpack = np.linspace(*thicks[i-1])
                angle_unpack = np.linspace(*angles[i])
                assert len(thick_unpack) == len(angle_unpack)
                for j in range(len(angle_unpack)):
                    self.epsO.append(self.materials[layers[i]][0])
                    self.epsE.append(self.materials[layers[i]][1])
                    self.mu.append(self.materials[layers[i]][2])
                    self.angle.append(angle_unpack[j])
                    self.thick.append(thick_unpack[j])
                    self.inter.append(self.inter[-1]+thick_unpack[j])
            else:
                raise NotImplementedError(
                    'Layers, thickness or angles samples mismatch.')
        self.epsO.append(self.materials[layers[-1]][0])
        self.epsE.append(self.materials[layers[-1]][1])
        self.mu.append(self.materials[layers[-1]][2])
        self.angle.append(angles[-1])

    def addLayer(self, *args):
        for arg in args:
            self.epsO.append(self.materials[arg][0])
            self.epsE.append(self.materials[arg][1])
            self.mu.append(self.materials[arg][2])

    def addThickness(self, *args):
        for arg in args:
            self.thick.append(arg)
            self.inter.append(self.inter[-1]+arg)

    def addAngle(self, *args):
        for arg in args:
            self.angle.append(arg)

    def mMatrix(self, i):
        M = np.zeros((self.N, 4, 4), dtype=np.complex)
        kzO = self.freqs*np.sqrt(self.epsO[i]*self.mu[i])*self.thick[i-1]
        kzE = self.freqs*np.sqrt(self.epsE[i]*self.mu[i])*self.thick[i-1]
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

    def aMatrix(self, i):
        a = np.zeros((self.N, 4, 4), dtype=np.complex)
        kzO = 1j*self.freqs*np.sqrt(self.epsO[i]*self.mu[i])*self.inter[i]
        kzE = 1j*self.freqs*np.sqrt(self.epsE[i]*self.mu[i])*self.inter[i]
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

    def ainvMatrix(self, i):
        ainv = np.zeros((self.N, 4, 4), dtype=np.complex)
        kzO = 1j*self.freqs*np.sqrt(self.epsO[i]*self.mu[i])*self.inter[i-1]
        kzE = 1j*self.freqs*np.sqrt(self.epsE[i]*self.mu[i])*self.inter[i-1]
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

    def tMatrix(self, i, f):
        T = np.matmul(self.aMatrix(i), R2(self.angle[i]))
        for k in range(i+1, f):
            T = np.matmul(R1(self.angle[k]-self.angle[k-1]), T)
            T = np.matmul(self.mMatrix(k), T)
        T = np.matmul(R1(self.angle[f]-self.angle[f-1]), T)
        T = np.matmul(self.ainvMatrix(f), T)
        T = np.matmul(R2(-self.angle[f]), T)
        return T

    def calculateCoeff(self, i=0, f=None):
        if f == None:
            f = len(self.mu) - 1
        T = self.getT(i, f)
        for efield in self.IA:
            IA = self.IA[efield]
            II = np.zeros((self.N), dtype=complex)  # incident intensity
            RI = np.zeros((self.N), dtype=complex)  # reflected intensity
            TI = np.zeros((self.N), dtype=complex)  # transmitted intensity
            A = np.zeros((self.N, 4, 4), dtype=np.complex)
            B = np.zeros((self.N, 4), dtype=np.complex)
            A[:, 0, 0] = T[:, 0, 1]
            A[:, 0, 1] = T[:, 0, 3]
            A[:, 0, 2] = -1
            A[:, 0, 3] = 0
            A[:, 1, 0] = T[:, 1, 1]
            A[:, 1, 1] = T[:, 1, 3]
            A[:, 1, 2] = 0
            A[:, 1, 3] = 0
            A[:, 2, 0] = T[:, 2, 1]
            A[:, 2, 1] = T[:, 2, 3]
            A[:, 2, 2] = 0
            A[:, 2, 3] = -1
            A[:, 3, 0] = T[:, 3, 1]
            A[:, 3, 1] = T[:, 3, 3]
            A[:, 3, 2] = 0
            A[:, 3, 3] = 0
            B[:, 0] = -T[:, 0, 0]*IA[:, 0] - T[:, 0, 2]*IA[:, 1]
            B[:, 1] = -T[:, 1, 0]*IA[:, 0] - T[:, 1, 2]*IA[:, 1]
            B[:, 2] = -T[:, 2, 0]*IA[:, 0] - T[:, 2, 2]*IA[:, 1]
            B[:, 3] = -T[:, 3, 0]*IA[:, 0] - T[:, 3, 2]*IA[:, 1]
            X = np.linalg.solve(A, B)
            RA = X[:, 0:2]
            TA = X[:, 2:4]
            for ii in range(self.N):
                II[ii] = np.vdot(IA[ii, :], IA[ii, :])
                RI[ii] = np.vdot(RA[ii, :], RA[ii, :])
                TI[ii] = np.vdot(TA[ii, :], TA[ii, :])
            self.coeff[f'{efield}({i},{f})'] = [
                IA, TA, RA, (TI/II).real, (RI/II).real, 1 - (TI/II).real - (RI/II).real]

    def TMM(self, i=0, f=None):
        if f == None:
            f = len(self.mu) - 1
        self.freqs = self.freqs * eV2THz  # freq from eV to THz
        self.T[f'T({i},{f})'] = self.tMatrix(i, f)

    def getT(self, i=0, f=None):
        if f == None:
            f = len(self.mu) - 1
        return self.T[f'T({i},{f})']

    def getCoeff(self, efield, i=0, f=None):
        if f == None:
            f = len(self.mu) - 1
        return self.coeff[f'{efield}({i},{f})']
