#!/bin/bash
#$ -l h_rt=06:00:00
#$ -cwd -V
#$ -N all_sims
#$ -j y
#$ -l np=16
#$ -l placement=optimal
# -t 1:5

source $HOME/.bashrc

#################################################################
################ Set the Parameters for the array ###############
#################################################################
COMMAND1='python run.py SAC --mpi' 
COMMAND2='python run.py gdf --mpi'

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
./configure.py compile sac --clean;

#### Run the CODE! ####
#time $COMMAND1
time $COMMAND2

###### I Done it now ########
rm -r $TMP_DIR
pushover -m "Job" ${exp_facs[i]} "Complete"
echo "Job" ${exp_facs[i]} "Complete"
