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
module load apps/python/conda
source activate mpi-sac

#Load MPI modules
module add mpi/gcc/openmpi/1.10.0

echo "SAC will run on the following nodes"
cat $PE_HOSTFILE

echo "Working DIR:"
echo $SGE_O_WORKDIR

echo "Run SAC:"
#time python $SGE_O_WORKDIR/run.py SAC --mpi

echo "Run GDF Translator:"
time python $SGE_O_WORKDIR/run.py gdf --mpi

echo "Job Complete"
