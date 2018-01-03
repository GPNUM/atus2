
[![Build Status](https://travis-ci.org/GPNUM/atus2.svg?branch=master)](https://travis-ci.org/GPNUM/atus2)

# **ATUS-2**
## **Introduction**
The ATUS-2 package is a collection of C++ programs designed to solve the non linear Schrödinger equation, especially the Gross--Pitaevskii equation. ATUS-2 is developed at ZARM (Center of Applied Space Technology and Microgravity, University of Bremen) and is supported by the German Space Agency (DLR) with funds provided by the Federal Ministry of Economics and Technology BMWi) due to an enactment of the German Bundestag under Grant No. DLR 50WM0942, 50WM1042, 50WM1342, 50WM1642.

The main purpose of ATUS-2 is to conduct numerical simulations of atom interferometers with Bose--Einstein condensates (BEC). At the moment the following light-matter interactions are supported for single as well as double species BECs (BEC-mixtures): Bragg, double Bragg and Raman. It solves the Cauchy initial value problem for the Gross--Pitaevskii equation.

The solvers are based on C++ template classes which are located in the include folder. These are instantiated for up to three spatial dimensions. Based on the form of the non linearity N the wave function can have an arbitrary number of internal degrees of freedom. All available programs are controlled through xml parameter files.

The code utilizes OpenMP and MPI for parallel computing.

## **Requirements and Dependencies**  

* make, cmake
* gcc, g++, gfortran
* Modules (recommended)
* doxygen (optional)
* Steel Bank Common Lisp (recommended for install script)
* gnuplot
* Paraview (optional)

The following packages are required and will be installed via the install script:

* openmpi 3.0.0
* gnu gsl 2.4
* muparser 2.2.5
* fftw 3.3.7
* lis 2.0.4
* vtk 8.0.1
* hdf5 1.10.1

The following third party packages are included in the source tree:

* pugixml (http://pugixml.org/)
* String Toolkit (http://www.partow.net/programming/strtk/index.html)
* CXXOPTS (https://github.com/jarro2783/cxxopts)


## **Installation**

It is highly recommended to use our install script which downloads all required packages and installs everything in $HOME/local/ by default. In order to run this script Steel Bank Common Lisp (sbcl) is required.

Further we recommended the use of environmental modules (http://modules.sourceforge.net/) for setting up all paths. Modules should be available via your Linux distribution. The install script also generates module files which are located in $HOME/local/modules. The search paths for the module files needs to be extended by adding your path to $MODULEPATH of your shell.

The binaries are installed in $HOME/bin. Make sure that this folder is added to $PATH of your shell.

## **Example 1 - double slit experiment**
Change to the sub folder xml/double_slit.

### Step 1: Preparation of the initial wave function
gen_psi_0 runme.xml

### Step 2: Real time propagation of the wave function
rt_solver runme.xml

### Step 3: View the result
Run the gnuplot script via: gnuplot plotall.gnuplot.

This script generates jpeg files of all stored intermediate time steps.

## **Example 2 - Bragg pulse with Rubidium BEC**
Change to the sub folder xml/bragg.

### Step 1: Preparation of the initial wave function
Run gen_psi_0 gauss.xml; this generates the initial wave function for the first internal state.

Run gen_psi_0 zero.xml; this generates the initial wave function for the second internal state.

### Step 2: Real time propagation of the wave function
Run bragg bragg_ad.xml

This command propagates the initial wave functions. It consists of two phases: The first is with Laser-matter interaction enabled; the second is free propagation. This generates binary output files.

### Step 3: View the result
Now the result can inspected with gnuplot.

The following gnuplot command

plot "<gpo3 7100.000_1.bin" w l

plots the final result of the last time step using the pipe program gpo3. By default the density is plotted.

The file Rabi_1_0.txt contains the Rabi oscillations. This can also be viewed with gnuplot.

## **Example 3 - Groundstate**
Change to the sub folder xml/groundstate.

### Step 1: Breed the ground state in a harmonic trap
Run sobmin groundstate.xml.

### Step 2: View the result
Use the following gnuplot commands

set pm3d map
splot "<gpo3 final.bin"

## **Example 4 - Groundstate of two coupled Gross--Pitaevskii equations**
Change to the sub folder xml/groundstate.

### Step 1: Breed the ground state in a harmonic trap
Run sobmin_2 groundstate_2.xml.

### Step 2: View the result
Use the following gnuplot command

p "<gpo3 final_1.bin" w l, "<gpo3 final_2.bin" w l
