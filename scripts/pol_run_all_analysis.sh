#!/bin/bash
#$ -l h_rt=15:00:00
#$ -cwd -V
#$ -N all_analysis
#$ -j y
#$ -l np=16
#$ -l placement=scatter
# -t 1:5

source $HOME/.bashrc

#################################################################
################ Set the Parameters for the array ###############
#################################################################

exp_facs=( 0.015 0.05 0.15 0.45 1.5 )

#################################################################
####################### Run the Script ##########################
#################################################################

#### Setup and Configure ####
#i=$((SGE_TASK_ID - 1))
i=2

BASE_DIR=$HOME/BitBucket/expansion-factor-paper/
TMP_DIR=$(mktemp -d --tmpdir=/nobackup/shesm/temp_run/)

cp -r $BASE_DIR $TMP_DIR
cd $TMP_DIR/expansion-factor-paper/

echo "exp" ${exp_facs[i]}
./configure.py set driver --exp_fac=${exp_facs[i]};
./configure.py print;
#./configure.py compile sac --clean;

#### Run the CODE! ####
tube_radii=( 'r10' 'r30' 'r60' )
for t in {1..2}
do
    echo ${tube_radfii[$t]}
    time python run.py analysis --mpi --tube-r=${tube_radii[$t]}
done

###### I Done it now ########
rm -r $TMP_DIR
pushover -m '"Job" ${exp_facs[i]} "Complete"'
echo "Job" ${exp_facs[i]} "Complete"
