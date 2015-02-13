#!/bin/bash
#$ -l h_rt=99:00:00
#$ -cwd
#$ -l arch=intel*
#$ -l mem=6G
#$ -pe openmpi-ib 16
#$ -N S3D_Slog
#$ -P mhd
#$ -q mhd.q
#$ -j y

#Set the Python virtualenv
source ~/.bashrc
workon vtk_hdf 

#Load MPI modules
module add mpi/pgi/openmpi/1.6.4

echo "SAC will run on the following nodes"
cat $PE_HOSTFILE

echo "Run SAC:"
time python /home/smq11sjm/period-paper/run.py SAC --mpi

echo "Run GDF Translator:"
time python /home/smq11sjm/period-paper/run.py gdf --mpi

echo "Job Complete"
/home/smq11sjm/.local/bin/pushover "Job Complete S3D_Slog"
