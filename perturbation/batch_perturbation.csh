#!/bin/csh

#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --error=batch_perturbation.%j.err
#SBATCH --output=batch_perturbation.%j.log
#SBATCH --exclusive

./batch_in_silico_perturbation.py
