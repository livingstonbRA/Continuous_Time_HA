#!/bin/bash
#SBATCH --job-name=combine2a
#SBATCH --output=/home/livingstonb/GitHub/Continuous_Time_HA/output/combine.out
#SBATCH --error=/home/livingstonb/GitHub/Continuous_Time_HA/output/combine.err
#SBATCH --partition=broadwl
#SBATCH --array=1
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8000

module load matlab/2019b

matlab -nodisplay < /home/livingstonb/GitHub/Continuous_Time_HA/code/batch/combine_runs.m
