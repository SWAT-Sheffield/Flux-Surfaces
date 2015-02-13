#!/bin/bash
#qsub mpi_run_mhd.sh;
#sleep 5s;
qsub mpi_analysis_r10.sh;
qsub mpi_analysis_r30.sh;
qsub mpi_analysis_r60.sh;
