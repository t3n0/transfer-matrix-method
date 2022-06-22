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
import matplotlib.pyplot as plt
from modules.material import Material
from modules.utils import input

data = input()

default = [0.0, 10.0, 1000]
materials = {'air': 1.0,
             'cnt': 'dielectrics/cnt',
             'sapp': [3.07, 0.07]}
efields = {"x_pol": 1.0,
           "y_pol": [1.0, "y"],
           "left_pol": [1.0, "left"],
           "THz_pulse": "fields/thz_pulse/"}

geometries = {
    "setup1": {
        "layers": ["air", "cnt", "air"],
        "angles": [0.0, [0.0, np.pi/2, 20], 0.0],
        "thicks": [[5.0, 5.0, 20]]
    }}

mat = Material(materials, efields, default)

mat.setGeometry(**geometries["setup1"])

mat.TMM()

mat.calculateCoeff()

ia, ra, ta, ti, ri, ai1 = mat.getCoeff('x_pol')
ia, ra, ta, ti, ri, ai2 = mat.getCoeff('y_pol')

plt.plot(mat.freqs, ai1)
plt.plot(mat.freqs, ai2)
plt.show()

# N = 20
# mat.addEpsE(1.0)
# mat.addEpsO(1.0)
# mat.addMu(1.0)
# mat.addAngle(0.0)
# for i in range(N):
#     mat.addEpsE(epsEpath)
#     mat.addEpsO(epsOpath)
#     mat.addMu(1.0)
#     mat.addAngle(+1*1*np.pi/2)
#     mat.addLayer(5.0)
# mat.addEpsE(1.0)
# mat.addEpsO(1.0)
# mat.addMu(1.0)
# mat.addAngle(0.0)

# # jones vectors
# E0 = 1.0
# x_pol     = E0*np.array([1, 0], dtype= complex)
# y_pol     = E0*np.array([0, 1], dtype= complex)
# left_pol  = E0*np.array([1, -1j], dtype= complex)
# right_pol = E0*np.array([1, 1j], dtype= complex)

#mat.addIncidentLight([1,2j], 'left')

# mat.TMM()

# i_pol1, t_pol1, r_pol1, T1, R1, A1 = mat.coeff(y_pol)
#t_pol2, r_pol2, T2, R2, A2 = mat.coeff(left_pol)


# plt.plot(mat.freqs/eV2THz, 0.0*np.ones((mat.N)), 'k--')
# plt.plot(mat.freqs/eV2THz, abs(t_pol1[:,1]), label='|Ey|')
# plt.plot(mat.freqs/eV2THz, abs(t_pol1[:,0]), label='|Ex|')
#plt.plot(mat.freqs/module.eV2THz, np.angle(t_pol1[:,1]) - np.angle(t_pol1[:,0]), label='ellipticity')
# plt.plot(mat.freqs/eV2THz, angle_between(t_pol1[:,1], t_pol1[:,0]), label='ellipticity2')
#plt.plot(mat.freqs/module.eV2THz, np.angle(t_pol1[:,0]), label='x')
#plt.plot(mat.freqs/module.eV2THz, np.angle(t_pol1[:,1]), label='y')
#plt.ylim(-np.pi/2, np.pi/2)
# plt.xlabel('Frequency (eV)')
# plt.ylabel('Transmitted fields')

#plt.plot(mat.freqs/module.eV2THz, A1, label=f'x pol. {jj}')
#plt.plot(mat.freqs/module.eV2THz, A2, label=f'y pol. {jj}')
#plt.plot(mat.freqs/module.eV2THz, A1-A2, label=f'{jj}')
#plt.ylim(-0.001, 0.009)
#plt.ylim(-0.02, 0.13)
#plt.ylim(-0.05, 0.6)
# plt.xlabel('Frequency (eV)')
# plt.ylabel('Absorption coeff')

# plt.plot(mat.freqs, A2 - A1, label='A')
# plt.plot(mat.freqs, T2 - T1, label='T')
# plt.plot(mat.freqs, R2 - R1, label='R')
# plt.legend()

# plt.show()
