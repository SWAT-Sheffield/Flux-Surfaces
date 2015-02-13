#!/bin/bash
#$ -l h_rt=05:00:00
#$ -l placement=scatter
#$ -cwd -V
#$ -N Slog_anal
#$ -l np=321
#$ -j y
#$ -t 1-3

#Target Tube Radii
radii=("r10" "r30" "r60")

#Set the Python virtualenv
source ~/.bashrc

echo "Run analysis script:"
~/BitBucket/period-paper/run.py analysis --tube-r=${radii[${SGE_TASK_ID}]} --mpi

echo "Job Complete"
pushover -m "Job '"$JOB_NAME"' - '${radii[${SGE_TASK_ID}]}' complete"
