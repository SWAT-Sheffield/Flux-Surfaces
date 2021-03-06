#!/bin/bash
#$ -l h_rt=15:00:00
#$ -cwd
#$ -l arch=intel*
#$ -l mem=6G
#$ -pe openmpi-ib 16
#$ -P mhd
#$ -q mhd.q
#$ -N all_periods
#$ -j y
#$ -t 1 

source $HOME/.bashrc
module load mpi/intel/openmpi/1.8.3
module load apps/python2-virtual/2.7.6
workon vtk6
export PATH=:/home/smq11sjm/.virtualenvs/vtk6/bin:$PATH
export LD_LIBRARY_PATH=/home/smq11sjm/.virtualenvs/vtk_hdf/lib/vtk-6.1/:$LD_LIBRARY_PATH
echo $(which python)
#################################################################
################ Set the Parameters for the array ###############
#################################################################
drivers=('Slog')

#################################################################
####################### Run the Script ##########################
#################################################################

#### Setup and Configure ####
i=$((SGE_TASK_ID - 1))

BASE_DIR=$HOME/GitHub/SWAT/Flux-Surfaces
TMP_DIR=$(mktemp -d --tmpdir=/fastdata/smq11sjm/temp_run/)

cp -r $BASE_DIR $TMP_DIR
cd $TMP_DIR/Flux-Surfaces/
pwd

echo ${periods[i]} ${amps[i]} "${fortamps[i]}"
#./configure.py set driver --period=${periods[i]} --amp=${amps[i]} --fort_amp="${fortamps[i]}";
./configure.py set driver --driver=${drivers[i]};
./configure.py print;
./configure.py compile sac --clean;

#### Run SAC ####
python ./run.py SAC --mpi

#### Run the CODE! ####
tube_radii=( 'r60' 'r30' 'r10' )
for tuber in "${tube_radii[@]}"
do
    echo $tuber
    python ./run.py analysis --mpi --tube-r=$tuber
done
###### I Done it now ########
rm -r $TMP_DIR
pushover -m "Job "${drivers[i]}" Complete"
echo "Job "${drivers[i]}" Complete"
