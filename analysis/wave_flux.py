#!/usr/bin/env python2
"""
This script re-calculates the Wave Energy Flux from the saved gdf files and 
outputs new flux gdf files and also saves new WaveFlux vtk files and re-calculates
the time distance arrays using the indices saved from the original run.

Usage:
    surface_analysis_mpi.py
"""
import os
import sys
import glob

import h5py
import numpy as np
import astropy.units as u
import yt.mods as ytm

import pysac.analysis
import pysac.io
import pysac.io.gdf_writer
import pysac.io.yt_fields

sys.path.append('../')
from scripts import sacconfig
cfg = sacconfig.SACConfig()

try:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    mpi = True
    mpi = mpi and (size != 1)
except ImportError:
    mpi = False
    rank = 0
#==============================================================================
# Config
#==============================================================================
cfg = sacconfig.SACConfig()
#Set defaults for no arguments
driver = cfg.driver
post_amp = cfg.amp
period = cfg.str_period

#Add the '_' to exp_fac
if cfg.exp_fac:
    exp_fac = '_' + cfg.str_exp_fac
else:
    exp_fac=''

#Make the pull paths
data_dir = cfg.data_dir
identifier = cfg.get_identifier()

gdf_path = cfg.gdf_dir

gdf_files = glob.glob(os.path.join(gdf_path, identifier+'_0*.gdf'))
gdf_files.sort()
timeseries = ytm.load(gdf_files)

if rank == 0:
    print "Configuration:"
    print 'driver:', driver
    print 'post_amp:', post_amp
    print 'period:', period
    print 'exp_fac:', exp_fac
    print 'data_dir:', data_dir
    print 'gdf_dir:', gdf_path

if rank == 0: #Prevents race condition where one processes creates the dir
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)

#Define a var to limit iterations, no limt = len(ts)
max_n = len(timeseries)

def path_join(filename):
    return os.path.join(data_dir,filename)

def glob_files(tube_r, search):
    files = glob.glob(os.path.join(cfg.data_dir,tube_r,search))
    files.sort()
    return files

#Divide up the time steps to each process
all_indices = np.arange(0, len(timeseries), 1, dtype=int)
if mpi:
    chunks = np.array_split(all_indices, size)
    rank_indices = comm.scatter(chunks, root=0)
else:
    rank_indices = all_indices

#Create line output arrays
#Number of indicies for this process:
number_ind = len(rank_indices)

def do_flux_step(k,i):
    print i
    ds = timeseries[i]
    
    f_wave = pysac.analysis.get_wave_flux_yt(ds)
    
    #Save waveflux to the gdf file:
    fname, ext = os.path.splitext(ds.parameter_filename)
    fprefix = fname[:-6]
    fnumber = fname[-6:]
    new_name = fprefix + '_fwave' + fnumber + ext

    f = h5py.File(new_name, mode='w')
    print i, new_name
    f = pysac.io.gdf_writer.create_file(f,
                                        {'ndim':3, 'nx':ds.domain_dimensions,
                                         't': ds.current_time},
                                         domain_left_edge=ds.domain_left_edge,
                                         domain_right_edge=ds.domain_right_edge)

    fwave_x = u.Quantity(f_wave[0,:,:,:], unit='W/m2')
    pysac.io.gdf_writer.write_field_u(f, fwave_x, 'wave_flux_x', 'x Compoment of Wave Energy Flux')
    fwave_y = u.Quantity(f_wave[1,:,:,:], unit='W/m2')
    pysac.io.gdf_writer.write_field_u(f, fwave_y, 'wave_flux_y', 'y Compoment of Wave Energy Flux')
    fwave_z = u.Quantity(f_wave[2,:,:,:], unit='W/m2')
    pysac.io.gdf_writer.write_field_u(f, fwave_z, 'wave_flux_z', 'z Compoment of Wave Energy Flux')
    f.close()

for k,i in enumerate(rank_indices):
    do_flux_step(k,i)