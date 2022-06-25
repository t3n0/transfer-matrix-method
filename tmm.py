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

from fileinput import filename
from modules.material import Material

if __name__ == '__main__':
    import os
    from modules.utils import input

    data = input()
    geometries = data['geometries']
    mat = Material(data['materials'], data['efields'], data['freq'])

    WORK_DIR = os.getcwd()
    MATS_DIR = os.path.join(WORK_DIR, data['flag'], 'materials')
    GEOS_DIR = []
    for geo in geometries:
        GEOS_DIR.append(os.path.join(WORK_DIR, data['flag'], geo))

    if not os.path.exists(MATS_DIR):
        os.makedirs(MATS_DIR)
    for geo_dir in GEOS_DIR:
        if not os.path.exists(geo_dir):
            os.makedirs(geo_dir)

    # plotting materials and saving to file
    for name in mat.materials:
        filename = os.path.join(MATS_DIR, name)
        mat.plotMaterial(name, filename)

    for ig, geo in enumerate(geometries):
        mat.setGeometry(**geometries[geo])
        mat.TMM()
        mat.calculateCoeff()

        # loop oover the efields?
        x_pol = mat.getCoeff('x_pol')
        ia, ta, ra, ti, ri, ai2 = mat.getCoeff('left_pol')

        # plotting geometry setup
        filename = os.path.join(GEOS_DIR[ig], geo)
        mat.plotGeometry(geo, filename)
        # plotting results and save to file
        #mat.plotResults()
        #mat.saveResults()
