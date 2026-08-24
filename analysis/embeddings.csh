#!/bin/csh

#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --error=embeddings.%j.err
#SBATCH --output=embeddings.%j.log
#SBATCH --exclusive

./extract_and_plot_cell_embeddings.py
