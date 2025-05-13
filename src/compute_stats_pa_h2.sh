#!/bin/sh
#=============================================#
#  script for launching compute_stats_pa_h2.R #
#=============================================#
#SBATCH -A wiser
#SBATCH --partition fast
#SBATCH --mem 90GB
#SBATCH --cpus-per-task 12
source /etc/profile.d/modules.sh
module load r/4.4.1
Rscript compute_stats_pa_h2.R
