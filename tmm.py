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

from modules.material import Material
import matplotlib.pyplot as plt

if __name__ == '__main__':
    import os
    from modules.utils import input

    data = input()
    geometries = data['geometries']
    mat = Material(data['materials'], data['efields'], data['freq'])

    WORK_DIR = os.getcwd()
    MATS_DIR = os.path.join(WORK_DIR, data['flag'], 'materials')
    GEOS_DIR = [os.path.join(WORK_DIR, data['flag'], geo) for geo in geometries]
    print(f'Current working directory:\n\t{WORK_DIR}')

    if not os.path.exists(MATS_DIR):
        os.makedirs(MATS_DIR)
    for geo_dir in GEOS_DIR:
        if not os.path.exists(geo_dir):
            os.makedirs(geo_dir)
    print(f'Output directories:\n\t{MATS_DIR}', end='')
    for i in range(len(GEOS_DIR)): print(f'\n\t{GEOS_DIR[i]}', end='')
    print()

    # plotting dielectric functions of all materials and saving to file
    for name in mat.materials:
        print(f'Saving dielectric functions to file: {name}')
        filename = os.path.join(MATS_DIR, name)
        mat.plotMaterial(name, filename)

    # loop over all geometry setups
    print(f'\n-------- Calculation starts --------\n')
    for ii, geo in enumerate(geometries):
        print(f'Set geometry: {geo}')
        mat.setGeometry(**geometries[geo])
        # plotting geometry setup
        filename = os.path.join(GEOS_DIR[ii], geo)
        mat.plotGeometry(geo, filename)

        print(f'\tCompute TMM matrix...')
        mat.TMM()
        print(f'\tCalculate coefficients...')
        mat.calculateCoeff()


        for ef in data['efields']:
            print(f'\t\tIncident fields: {ef}')
            coeff = mat.getCoeff(ef)
            # plotting results and save to file
            filename = os.path.join(GEOS_DIR[ii], ef)
            mat.plotResults(ef, filename)
