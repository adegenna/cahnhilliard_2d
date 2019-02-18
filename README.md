# Introduction
Solver for the 2D Cahn Hilliard equation. Zero flux boundary conditions, square domain.

# Solver
First, edit the filepaths in src/input_driver.dat to reflect what is on your machine. Then, run:

python driver.py input_driver.dat

# Visualization
First, edit the filepaths in the visualization/plotContourState.py script to reflect what is on your machine. Then, run:

python plotContourState.py

# Example Output
Here are some example state snapshots from a simulation run with epsilon = 0.1, dt = 5e-5
<img src="https://github.com/adegenna/cahnhilliard_2d/blob/master/ch2d.png">
