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

geometries = data['geometries']

mat = Material(data['materials'], data['efields'], data['freq'])

mat.setGeometry(**geometries["setup1"])
mat.TMM()
mat.calculateCoeff()

ia, ra, ta, ti, ri, ai1 = mat.getCoeff('x_pol')
ia, ra, ta, ti, ri, ai2 = mat.getCoeff('left_pol')

plt.plot(mat.freqs, ai1)
plt.plot(mat.freqs, ai2)

mat.setGeometry(**geometries["setup2"])
mat.TMM()
mat.calculateCoeff()

ia, ra, ta, ti, ri, ai1 = mat.getCoeff('x_pol')
ia, ra, ta, ti, ri, ai2 = mat.getCoeff('left_pol')

plt.plot(mat.freqs, ai1)
plt.plot(mat.freqs, ai2)
plt.show()
