# Introduction
Solver for the 2D Cahn Hilliard equation:

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/thermal/cheqn.gif">

Domain: 2D rectangular

Boundary Conditions: periodic, Neumann, Dirichlet, or mixed

Some of the coefficients of this equation are functions of temperature:

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/thermal/eps2_thermal.gif">

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/thermal/sigma_thermal.gif">

Temperature itself can be evolved as a spatial field through a thermal diffusion equation that is one-way coupled to the Cahn-Hilliard dynamics:

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/thermal/thermal_eqn.gif">

# cpp/
This subdirectory contains all C++ source (in `src/`) as well as Python Swig wrappers (in `swig`).

## Building: Pure C++ 
Pure C++ executables may be built with the included `CMakeLists.txt` as follows:

```shell
cd [/PATH/TO]/cahnhilliard_2d/cpp
mkdir build/
cd build/
cmake ../
make
```

The resulting executable may be run:
```shell
./ch2d
```

## Building: C++/Swig
A separate Makefile is included to wrap the C++ source code into a Python module that can be imported and used with Python. To make, do:

```shell
cd [/PATH/TO]/cahnhilliard_2d/cpp
make
```

This will compile the C++ code into a .so shared object and move it into the `swig/` subdirectory. The `swig/` subdirectory provides an example driver demonstrating how to use this:

```shell
cd [/PATH/TO]/cahnhilliard_2d/cpp/swig
python driver.py
```

# python/
This subdirectory contains earlier prototyping code in pure Python, but support for it has been discontinued in favor of the C++/Swig solution.

# Visualization
Some lightweight Python scripts are provided in `visualization/` for convenience. You will have to edit the necessary options/filepaths to be consistent with the state files you are trying to read/display.

# Example Output
Here are some example state snapshots from 4 different simulations. The first three show the effect of increasing the biharmonic coefficient. The last snapshot is taken from the same system as that which produced the second snapshot, but with more noise added to the dynamics.

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/ch2d.png">

These are two steady-states achieved with differing values of `m` (the left has `m = 0.5`; the right has `m = 0`):

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/thermal/ch_nonthermal.png" width="200" height="200"> <img src="https://github.com/adegenna/cahnhilliard_2d/blob/thermal/ch_nonthermal_2.png" width="200" height="200">

This is a temperature-dependent simulation with thermal diffusion present (the top is the concentration field; the bottom is the temperature field):

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/thermal/ch_thermal.png" width="200" height="400">
