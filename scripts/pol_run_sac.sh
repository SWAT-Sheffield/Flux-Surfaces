#!/bin/bash
#$ -l h_rt=10:00:00
#$ -cwd -V
#$ -N Slog_period
#$ -l np=16
#$ -j y

period=270.0
amp=A10
fort_amp=10.d0

#Set the Python env
source $HOME/.bashrc

#Make sure we are in the right dir
cd $HOME/BitBucket/period-paper/

#configure the repo
python configure.py set driver --period=$period --amp=$amp --fort_amp=$fort_amp
python configure.py compile SAC --clean

echo "SAC will run on the following nodes"
#cat $PE_HOSTFILE

echo "Run SAC:"
#time python run.py SAC --mpi

echo "Run GDF Translator:"
#time python run.py gdf --mpi

echo "Job Complete"
pushover -m "Job '"$JOB_NAME"' Complete"
