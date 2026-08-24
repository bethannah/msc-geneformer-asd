#!/bin/csh

#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --error=tokenizing.%j.err
#SBATCH --output=tokenizing/%j.log
#SBATCH --exclusive

./tokenization.py
