#!/bin/sh
#=============================================#
#  script for launching compute_stats_pa_h2.R #
#=============================================#
### Requirements
#SBATCH --partition=p01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=90
#SBATCH --cpus-per-task=8
Rscript compute_stats_pa_h2.R
