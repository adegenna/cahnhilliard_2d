# Introduction

## Concentration dynamics
The main field that evolves is relative concentration of the A phase at the mesoscale of a block copolymer, which occurs through a modified 2D Cahn Hilliard equation:

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/cheqn.gif">

Domain: 2D rectangular

Boundary Conditions: periodic, Neumann, Dirichlet, or mixed

## Thermal dynamics

Some of the coefficients of this equation may be chosen to be a function of temperature, based on scaling arguments:

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/eps2_thermal.gif">

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/sigma_thermal.gif">

Temperature itself can be evolved as a spatial field through a thermal diffusion equation that is one-way coupled to the Cahn-Hilliard dynamics:

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/thermal_eqn.gif">

## Numerical method

Spatial discretization is uniform finite difference. Temporal evolution is handled by the `boost::numeric::odeint` library. A variety of explicit marching options are possible (e.g., RK4, with or without adaptive steps). The explicit timestep is set with reference to the biharmonic timescale of the problem (which should be the stiffest linear timescale present).

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

# Parallel Scaling
The solver is parallelized using shared memory with OpenMP. Here are some small scaling studies, using no thermal behavior and constant CH coefficients. In what follows:

* `nx` = spatial resolution
* `epsilon` = value of biharmonic coefficient
* `n_dtbiharm` = temporal length of total simulation, referenced to the length of the biharmonic timescale
* `n_core` = number of CPU cores
* `time (sec)` = total execution time, in seconds

## BloodMeridian (6 cores, Intel(R) Core(TM) i7-6800K CPU @ 3.40GHz)

### Pure C++
| `nx`          | `epsilon`     | `n_dtbiharm`  | `n_core`      | `time (sec)`  |
| ------------- |:-------------:| -------------:| ------------- |:-------------:|
| 128           | 0.01          | 300           | 1             | 305           |
| 128           | 0.01          | 300           | 3             | 41            |
| 128           | 0.01          | 300           | 6             | 36            |

### C++/Swig
| `nx`          | `epsilon`     | `n_dtbiharm`  | `n_core`      | `time (sec)`  |
| ------------- |:-------------:| -------------:| ------------- |:-------------:|
| 128           | 0.01          | 300           | 1             | 255           |
| 128           | 0.01          | 300           | 3             | 54            |
| 128           | 0.01          | 300           | 6             | 51            |

## HPC1 (8 cores, Intel(R) Xeon(R) CPU E5-2670 0 @ 2.60GHz)

### Pure C++
| `nx`          | `epsilon`     | `n_dtbiharm`  | `n_core`      | `time (sec)`  |
| ------------- |:-------------:| -------------:| ------------- |:-------------:|
| 128           | 0.01          | 300           | 1             | 309           |
| 128           | 0.01          | 300           | 4             | 52            |
| 128           | 0.01          | 300           | 8             | 34            |

# Example Output
Here are some example state snapshots from 4 different simulations. The first three show the effect of increasing the biharmonic coefficient. The last snapshot is taken from the same system as that which produced the second snapshot, but with more noise added to the dynamics.

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/ch2d.png">

These are two steady-states achieved with differing values of `m` (the left has `m = 0.5`; the right has `m = 0`):

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/ch_nonthermal.png" width="200" height="200"> <img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/ch_nonthermal_2.png" width="200" height="200">

This is a temperature-dependent simulation with thermal diffusion present (the top is the concentration field; the bottom is the temperature field):

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/ch_thermal.png" width="200" height="400">
