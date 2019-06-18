# Introduction
Solver for the 2D Cahn Hilliard equation:

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/cheqn.gif">

Zero flux boundary conditions, square domain.

# Solver
First, edit the filepaths in src/input_driver.dat to reflect what is on your machine. Then, run:

python driver.py input_driver.dat

# Visualization
First, edit the filepaths in the visualization/plotContourState.py script to reflect what is on your machine. Then, run:

python plotContourState.py

# Example Output
Here are some example state snapshots from 4 different simulations. The first three show the effect of increasing the biharmonic coefficient. The last snapshot is taken from the same system as that which produced the second snapshot, but with more noise added to the dynamics.

<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/ch2d.png">
