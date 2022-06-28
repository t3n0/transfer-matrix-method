# Transfer matrix method tool (tmm-tool)

[![GitHub Release Date](https://img.shields.io/github/release-date/t3n0/transfer-matrix-method)](https://github.com/t3n0/transfer-matrix-method/releases/latest)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/t3n0/transfer-matrix-method)](https://github.com/t3n0/transfer-matrix-method/releases/latest)
[![GitHub all releases](https://img.shields.io/github/downloads/t3n0/transfer-matrix-method/total)](https://github.com/t3n0/transfer-matrix-method/releases/latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

![tmm-tool](./output_sample/setup2/setup2.png)

This tool implements the transfer matrix method (TMM) to describe the propagation of electromagnetic (EM) waves in stacked layered materials at normal incidence. It uses a 4x4 matrix formalism which can account for the propagation of the EM fields with both x and y components. The EM waves are described using the [Jones vector formalism](https://en.wikipedia.org/wiki/Jones_calculus). This allows the user to study the behaviour of the system with respect to the polarisation of light. Also, any layer can be arbitrary rotated by a given angle, thus allowing to study complex responses such as linear and/or circular [dichroism](https://en.wikipedia.org/wiki/Dichroism).

## Installation

### Preferred installation: `pip`
The preferred installation method is by using `pip`. This will allow the package to be installed in Windows, Unix and MacOS. In the release page, simply download the [.zip file](https://github.com/t3n0/transfer-matrix-method/releases/latest) and extract it at your favorite location. Then move into the folder and run

`pip install .`

This will install two things:
 - the `tmm-tool` python package
 - and the `tmm-tool` shell command.
 
If you are using anaconda as you python environment, you can check the succesfull installation of the package by typing `conda list`.

To uninstall type `pip uninstall tmm-tool`.

### Stand-alone installation
If you are not confortable in using `pip` and prefer not to install the system-wide `tmm-tool` shell command, simply copy-paste the source and run the tool from within the base folder.

### Requirements
This tool requires:
- `python >=3.5`
- `numpy`
- `matplotlib`

## Usage

There are two ways of using `tmm-tool`:
- as a shell command via `~$ tmm-tool <my_input_file>`
- as a python package via `from tmm.material import Material`

### Shell command usage
Invoking the `tmm-tool` from the command line is the easiest way to start using the tool. Simply log into the terminal and type

`tmm-tool <my_input_file>` (for the pip installation case), or

`python tmm-tool <my_input_file>` (stand-alone installation, from the base directory).

The tool requires a file `<my_input_file>` in `json` format. An example of this is given in [input.txt](./input.txt) and a thorough description of the `json` (key, value) pairs is given in [here](not-yet).

The software will (in order):
- initialise the optical constants of the selected materials
- initialise the amplitudes of the incident electric fields
- construct the heterostructures (layers) given in input
- compute the total reflection, transmission, absorption coefficients (and more)
- output the data to the destination folders

At the moment, the command line code only computes the **total** transmission, reflection or absorption of the system.
In other words, the code computes the transfer matrix $T = T_{if}$ between the initial (i) and final (f) layers.
If the user wants to compute the coefficients for the intermidiate layers (e.g. as in a *echo removal* calculation), the python package (see below) should be used. This feature will be implemented in the next versions.

### Python package usage

To use `tmm-tool` as a package, you simply need to import the `Material` class

`from tmm.material import Material`


## Output

## Roadmap

The project is stable enough to deserve a release and offers the basic capabilities of the transfer matrix method. However, several improvements are on the making and will include:
 - a more efficient and versatile input file for scripting;
 - a routine to perform the removal of echos between layers;
 - a routine to calculate linear and circular dichroism;
 - the possibility of using different unit of measures (wavelength vs. frequency vs. energy);
 - an integration with the [RefractiveIndex](https://refractiveindex.info/) database, which provides the necessary optical constants.

## Support

For any problems, questions or suggestions, please contact me at tenobaldi@gmail.com.

## Authors and acknowledgment

The development of this tool is proudly powered by [me](https://github.com/t3n0).
Also, please consider citing the relevant literature if you are going to use this tool:
 - [Phys. Rev. B 104, 155437 (2021)](https://doi.org/10.1103/PhysRevB.104.155437)
 - [J Infrared Milli Terahz Waves 42, 1142–1152 (2021)](https://doi.org/10.1007/s10762-021-00815-5)

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
