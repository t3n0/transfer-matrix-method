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

import matplotlib.pyplot as plt
import numpy as np
import os
from modules.utils import *


SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


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
        self.freqs = self.freqs * eV2THz  # freq from eV to THz

    def setEpsMu(self, epsmu):
        self.epsO.append(epsmu[0])
        self.epsE.append(epsmu[1])
        self.mu.append(epsmu[2])

    def initGeometry(self):
        self.colors = {}
        for i, mat in enumerate(self.materials):
            self.colors[mat] = mycolors(i, len(self.materials))
        self.epsO = []
        self.epsE = []
        self.mu = []
        self.thick = []
        self.inter = [0.0]
        self.angle = []
        self.labelPos = []

    def setGeometry(self, layers, thicks, angles):
        self.colorlist = []
        self.layers = layers
        assert len(layers) == (len(thicks) + 2) == len(angles)
        self.initGeometry()
        self.setEpsMu(self.materials[layers[0]])
        self.colorlist.append(self.colors[layers[0]])
        self.angle.append(angles[0])
        for i in range(1, len(layers)-1):
            if (isinstance(thicks[i-1], int) or isinstance(thicks[i-1], float)) and (isinstance(angles[i], int) or isinstance(angles[i], float)):
                self.setEpsMu(self.materials[layers[i]])
                self.angle.append(angles[i])
                self.thick.append(thicks[i-1])
                self.labelPos.append(self.inter[-1] + self.thick[i-1]/2)
                self.inter.append(self.inter[-1]+self.thick[i-1])
                self.colorlist.append(self.colors[layers[i]])
            elif isinstance(thicks[i-1], list) and isinstance(angles[i], list):
                thick_unpack = np.linspace(*thicks[i-1])
                angle_unpack = np.linspace(*angles[i])
                assert len(thick_unpack) == len(angle_unpack)
                self.labelPos.append(self.inter[-1] + sum(thick_unpack)/2)
                for j in range(len(angle_unpack)):
                    self.setEpsMu(self.materials[layers[i]])
                    self.angle.append(angle_unpack[j])
                    self.thick.append(thick_unpack[j])
                    self.inter.append(self.inter[-1]+thick_unpack[j])
                    self.colorlist.append(self.colors[layers[i]])
            else:
                raise NotImplementedError(
                    'Layers, thickness or angles samples mismatch.')
        self.setEpsMu(self.materials[layers[-1]])
        self.colorlist.append(self.colors[layers[-1]])
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
        self.T[f'T({i},{f})'] = self.tMatrix(i, f)

    def getT(self, i=0, f=None):
        if f == None:
            f = len(self.mu) - 1
        return self.T[f'T({i},{f})']

    def getCoeff(self, efield, i=0, f=None):
        if f == None:
            f = len(self.mu) - 1
        return self.coeff[f'{efield}({i},{f})']

    def plotGeometry(self, title='', path=None):
        fig = plt.figure(figsize=(8, 4), dpi=300)
        ax1 = fig.add_axes([0.1, 0.4, 0.85, 0.45])
        fig.suptitle(title)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_xlim(self.inter[0], self.inter[-1])
        ax1.set_ylim(-1, 1)
        for i in range(len(self.inter)-1):
            xs = [self.inter[i], self.inter[i],
                  self.inter[i+1], self.inter[i+1]]
            ys = [-1, 1, 1, -1]
            c = self.colorlist[i+1]
            change = self.angle[i+1]/(max(self.angle)+1e-5)
            ax1.fill(xs, ys, color=(c[0], c[1]*change, c[2], 0.7))
        pos = 0
        y1 = 1.25
        y = y2 = 1.1
        for i in range(1, len(self.layers)-1):
            if self.labelPos[i-1] - pos < (self.labelPos[-1] - self.labelPos[0])/10:
                y = y1
                y1 = y2
                y2 = y
            pos = self.labelPos[i-1]
            ax1.text(self.labelPos[i-1], y,
                     self.layers[i], va='center', ha='center')
        ax1.text(-0.01*self.inter[-1], 0.0, self.layers[0],
                 va='center', ha='right', rotation=90)
        ax1.text(1.01*self.inter[-1], 0.0, self.layers[-1],
                 va='center', ha='left', rotation=90)

        ax2 = fig.add_axes([0.1, 0.2, 0.85, 0.2])
        for i in range(len(self.inter)-1):
            ax2.hlines(self.angle[i+1], self.inter[i], self.inter[i+1])
        ax2.set_xlim(self.inter[0], self.inter[-1])
        ax2.set_ylim(-0.05*max(self.angle), 1.05*max(self.angle))
        ax2.set_xlabel('z (nm)')
        ax2.set_ylabel('Θ (rad)')
        if path == None:
            plt.show()
        else:
            fig.savefig(path)
            plt.close()

    def plotMaterial(self, name, path=None):
        fig = plt.figure(figsize=(8, 6), dpi=300)
        ax1 = fig.add_axes([0.1, 0.1, 0.85, 0.85])
        fig.suptitle(name)
        ax1.plot(self.freqs/eV2THz, self.materials[name][0].real/eps_vac)
        ax1.plot(self.freqs/eV2THz, self.materials[name][0].imag/eps_vac)
        ax1.plot(self.freqs/eV2THz, self.materials[name][1].real/eps_vac)
        ax1.plot(self.freqs/eV2THz, self.materials[name][1].imag/eps_vac)
        ax1.set_xlabel('Frequency (eV)')
        ax1.set_ylabel('Dielectric function')
        if path == None:
            plt.show()
        else:
            plt.savefig(path)
            header = 'freq(eV) epsO1 epsO2 epsE1 epsE2'
            save_file(self.freqs/eV2THz, self.materials[name][0]/eps_vac,
                      self.materials[name][1]/eps_vac, path=path+'.txt', header=header)
            plt.close()

    def plotResults(self, efield, i=0, f=None):
        if f == None:
            f = len(self.mu) - 1
        fig = plt.figure(figsize=(8, 6), dpi=300)
        ax1 = fig.add_axes([0.1, 0.1, 0.4, 0.2])
        ax2 = fig.add_axes([0.1, 0.3, 0.4, 0.2])
        ax3 = fig.add_axes([0.1, 0.5, 0.4, 0.2])
        ax4 = fig.add_axes([0.6, 0.1, 0.4, 0.2])
        ax5 = fig.add_axes([0.6, 0.3, 0.4, 0.2])
        ax6 = fig.add_axes([0.6, 0.5, 0.4, 0.2])
        fig.suptitle(efield)
        #[f'{efield}({i},{f})']
        plt.show()
