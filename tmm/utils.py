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
import sys
import json
from scipy.interpolate import interp1d
from matplotlib.cm import get_cmap

# some physical constants
cspeed = 299792.458       # nm/ps
eps_vac = 8.8541878128e-9  # pF/nm
mu_vac = 1.25663706212e-3  # ps^2/pF/nm
eV2THz = 241.79893        # 1eV = 241.79893 THz


def mycolors(i, n):
    x = np.linspace(0.1, 0.9, n)
    cmap = get_cmap('rainbow')
    rgba = cmap(x)
    return rgba[i]


def input():
    args = sys.argv
    if len(args) < 2:
        print('No input file given. Exit.')
        sys.exit()
    else:
        input_file = args[1]
        with open(input_file, 'r') as f:
            data = json.load(f)
        return data


def freqMesh(materials, efield, default=[0.0, 10.0, 1000]):
    WORK_DIR = os.getcwd()
    freqs = []
    freqs.append([default[0], default[1]])
    for mat in materials:
        if isinstance(materials[mat], str):
            freqO, _ = load_file(os.path.join(
                WORK_DIR, materials[mat], 'epsO.txt'))
            freqE, _ = load_file(os.path.join(
                WORK_DIR, materials[mat], 'epsE.txt'))
            freqs.append(freqO)
            freqs.append(freqE)
    for ef in efield:
        if isinstance(efield[ef], str):
            freqEx, _ = load_file(os.path.join(WORK_DIR, efield[ef], 'Ex.txt'))
            freqEy, _ = load_file(os.path.join(WORK_DIR, efield[ef], 'Ey.txt'))
            freqs.append(freqEx)
            freqs.append(freqEy)
    freqmin, freqmax = minmax(*freqs)
    return np.linspace(freqmin, freqmax, default[2])


def angle_between(vec1, vec2):
    # angle between two complex vectors, without log branch cut line
    # i.e. the angle is not restricted between [-pi, pi]
    # the resulting curve is continous
    angle1 = np.angle(vec1)
    angle2 = np.angle(vec2)
    for i in range(len(angle1)-1):
        angle1[i+1] = angle1[i+1] + \
            np.rint((angle1[i] - angle1[i+1])/2/np.pi)*2*np.pi
        angle2[i+1] = angle2[i+1] + \
            np.rint((angle2[i] - angle2[i+1])/2/np.pi)*2*np.pi
    res = angle1 - angle2
    return res


def R1(theta):
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
    R = np.zeros((4, 4))
    for i in range(0, 4):
        R[i, i] = np.cos(theta)
    R[0, 2] = -np.sin(theta)
    R[1, 3] = +np.sin(theta)
    R[2, 0] = +np.sin(theta)
    R[3, 1] = -np.sin(theta)
    return R


def R2(theta):
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
    R = np.zeros((4, 4))
    for i in range(0, 4):
        R[i, i] = np.cos(theta)
    R[0, 2] = -np.sin(theta)
    R[1, 3] = -np.sin(theta)
    R[2, 0] = +np.sin(theta)
    R[3, 1] = +np.sin(theta)
    return R


def save_file(*args, path, header=''):
    data = []
    for arg in args:
        if (arg.dtype == float) or (arg.dtype == int):
            data.append(arg)
        elif (arg.dtype == complex):
            data.append(arg.real)
            data.append(arg.imag)
    data = np.array(data)
    np.savetxt(path, data.T, header=header)

def load_file(path):
    data = np.loadtxt(path)
    x = data[:, 0]
    if data.shape[1] == 2:
        y = interp1d(x, data[:, 1])
    else:
        y = interp1d(x, data[:, 1] + 1j*data[:, 2])
    return x, y


def minmax(*args):
    mins = []
    maxs = []
    for arg in args:
        mins.append(min(arg))
        maxs.append(max(arg))
    return max(mins), min(maxs)


def fill(vec, freqs):
    for i in range(0, len(vec)):
        if callable(vec[i]):
            vec[i] = vec[i](freqs)
        else:
            vec[i] = vec[i]*np.ones(len(freqs), dtype=np.complex)
    return np.array(vec)
