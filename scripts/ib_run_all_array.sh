#!/bin/bash
#$ -l h_rt=15:00:00
#$ -cwd
#$ -l arch=intel*
#$ -l mem=6G
#$ -pe openmpi-ib 12
#$ -P mhd
# -q mhd.q
#$ -N all_periods
#$ -j y
#$ -t 3
#1-11

source $HOME/.bashrc
module load mpi/pgi/openmpi/1.6.4
module load apps/python2-virtual/2.7.6
workon vtk_hdf
export PATH=:/home/smq11sjm/.virtualenvs/vtk_hdf/bin:$PATH
export LD_LIBRARY_PATH=/home/smq11sjm/.virtualenvs/vtk_hdf/lib/vtk-5.10/:$LD_LIBRARY_PATH
echo $(which python)
#################################################################
################ Set the Parameters for the array ###############
#################################################################
periods=( 30.0 60.0 90.0 120.0 150.0 180.0 210.0 240.0 270.0 300.0 330.0 )
amps=( A20r2 A20 A20r2-3 A10r2 A4r10 A20-r3 A20r2-7 A10 A20-3r2 A4r5 A20r2-11)
fortamps=( "'20.d0 * SQRT(2.d0)'" "'20.d0'" "'20.d0 * SQRT(2.d0 / 3.d0)'"
"'10.d0 * SQRT(2.d0)'" "'4.d0 * SQRT(10.d0)'" "'20.d0 * SQRT(3.d0)'" "'20.d0 * SQRT(2.d0 / 7.d0)'"
"'10.d0'" "'20.d0 / 3.d0 * SQRT(2.d0)'" "'4.d0 * SQRT(5.d0)'" "'20.d0 * SQRT(2.d0 / 11.d0)'")

#################################################################
####################### Run the Script ##########################
#################################################################

#### Setup and Configure ####
i=$((SGE_TASK_ID - 1))

BASE_DIR=$HOME/BitBucket/period-paper/
TMP_DIR=$(mktemp -d --tmpdir=/fastdata/smq11sjm/temp_run/)

cp -r $BASE_DIR $TMP_DIR
cd $TMP_DIR/period-paper/
pwd

echo ${periods[i]} ${amps[i]} "${fortamps[i]}"
./configure.py set driver --period=${periods[i]} --amp=${amps[i]} --fort_amp="${fortamps[i]}";
./configure.py print;
./configure.py compile sac --clean;

#### Run the CODE! ####
tube_radii=( 'r60' 'r30' 'r10' )
for tuber in "${tube_radii[@]}"
do
    echo $tuber
    python ./run.py analysis --mpi --tube-r=$tuber
done
###### I Done it now ########
rm -r $TMP_DIR
pushover -m "Job "${periods[i]}" Complete"
echo "Job "${periods[i]}" Complete"
