#!/bin/bash
#
#SBATCH --job-name=ch2d_petsc
#SBATCH --output=ch2d_petsc.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=3:00:00
#SBATCH --mem=64gb
# SBATCH --gres=gpu:1
# SBATCH -C k40

pwd; hostname; date
module load anaconda3
module load openmpi
source activate aeolus
python3 /home/adegennaro/cahnhilliard_2d/cpp/implicit_petsc/examples/driver.py
source deactivate
