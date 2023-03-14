#! /bin/sh
#SBATCH --cpus-per-task=6
#SBATCH --output results/example_R-%j.out
#SBATCH --error results/example_R-%j.err
# this assumes that the results directory exists already

srun /usr/bin/env Rscript --no-restore --no-save example_R.R
wait
