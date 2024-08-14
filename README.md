# THMC_FVM_Phreeqc
MRST: A sequential algorithm for THMC coupled modeling  based on the unified FVM framework
==============================================================================================

`THMC_FVM_Phreeqc` is is a finite volume method (FVM)-based coupled model designed to simulate thermal (T), hydraulic (H), mechanical (M), and chemical (C) processes in fractured reservoirs.  The model is developed using the Matlab Reservoir Simulation Toolbox (MRST) provided by SINTEF ICT (http://www.sintef.no/projectweb/mrst/). It employs the finite volume method (FVM) to solve fluid flow, heat transfer, and chemical reactions, while mechanical deformation is handled using the extended finite volume method (XFVM).The chemical reaction module integrates PHREEQC functionalities to simulate multicomponent reactions under high-temperature and high-pressure conditions.

Currently available features include:

* 2D Embedded discrete fracture modeling (EDFM)
* THMC coupled modeling based on the unified FVM framework
* Multicomponent reactions under high-temperature and high-pressure conditions

## Installation and Usage

The code has been tested and developed using `MRST 2021b`. Newer versions of MRST may be compatible.

It requires [iphreeqc-3.7.3-15968.tar.gz](https://water.usgs.gov/water-resources/software/PHREEQC/iphreeqc-3.7.3-15968.tar.gz)  to be installed on a Linux system before using the PHREEQC interface, and it also needs the 'geothermal.dat' file.

### Installation  ###

1. Install [MRST2021b](http://www.mrst.no) on a Linux system (CentOS 7.9 is used in this example).
2. Install [iphreeqc-3.7.3-15968.tar.gz](https://water.usgs.gov/water-resources/software/PHREEQC/iphreeqc-3.7.3-15968.tar.gz) on the Linux system.
3. Provide the 'geothermal.dat' database file.
4. Add all folders from `THMC_FVM_Phreeqc` to the `examples` folder under the `modules` directory in MRST.

Once MRST and this module are installed, this module can be used like any other MRST module.

### Getting started ###

To run examples in `THMC_FVM_Phreeqc` :

1. Before running any script that relies on the repository, MRST must be started. This is done by running the file `startup.m`, which is located inside your MRST directory.
2. Add all folders from `THMC_FVM_Phreeqc` to the path by selecting "Add to Path" -> "Selected Folders and Subfolders."
3. Double-click and run the case script in `THMC_XFVM_Phreeqc.m`.

## Project Overview

This project extends the capabilities of the MATLAB Reservoir Simulation Toolbox (MRST), which was originally developed by TU Delft and SINTEF ICT. The original software is licensed under the GNU General Public License, and this extended version is distributed under the same license.

Please note that this software also incorporates additional functions and libraries, which are subject to their respective licenses. Detailed information regarding these licenses can be found in the comments within the code.

## License

This project is licensed under the GNU General Public License along with MRST—see the [LICENSE]([GNU General Public License (GPL)](http://www.gnu.org/licenses/gpl.html).) file for details. The program also incorporates functionalities related to chemical reactions (from PHREEQC) and mechanical deformation, which have significantly contributed to the project. These contributions are acknowledged and recognized, though they are not individually listed here. For details on the licenses of these and other incorporated programs, please refer to the corresponding function descriptions within the code.

## Disclaimer

This software is provided "as is", without any warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and noninfringement. In no event shall the authors be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the software or the use or other dealings in the software.

## Reference

* K.-A. Lie. [An Introduction to Reservoir Simulation Using MATLAB/GNU Octave: User Guide for the MATLAB Reservoir Simulation Toolbox (MRST)](https://www.cambridge.org/core/books/an-introduction-to-reservoir-simulation-using-matlabgnu-octave/F48C3D8C88A3F67E4D97D4E16970F894). Cambridge University Press, 2019. [www.cambridge.org/9781108492430](https://www.cambridge.org/9781108492430)
* K.-A. Lie, S. Krogstad, I. S. Ligaarden, J. R. Natvig, H. M. Nilsen, and B. Skaflestad. [Open source MATLAB implementation of consistent discretisations on complex grids](http://folk.ntnu.no/andreas/papers/mrst-comg.pdf). Comput. Geosci., Vol. 16, No. 2, pp. 297-322, 2012. DOI: [10.1007/s10596-011-9244-4](http://dx.doi.org/10.1007/s10596-011-9244-4)
* Parkhurst D L, Appelo C A J. Description of input and examples for PHREEQC version 3—a computer program for speciation, batch-reaction, one-dimensional transport, and inverse geochemical calculations[J]. US geological survey techniques and methods, 2013, 6(A43): 497.https://pubs.usgs.gov/tm/06/a43/
* Zhang GR, Lu P, Zhang YL, Tu K, Zhu C. 2020. SUPPHREEQC: A program for generating customized PHREEQC thermodynamic datasets from SUPCRTBL and extending calculations to elevated pressures and temperatures. Computers & Geosciences, 143: 104560.  https://hydrogeochem.earth.indiana.edu 
* Deb R. Numerical modeling of fluid injection induced shear failure, tensile opening and flow-mechanics coupling[D]. ETH Zurich, 2018. https://doi.org/10.3929/ethz-b-000303929


