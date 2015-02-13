#! /usr/bin/env python2
"""
This script re-creates the vtk files with the surface fluxes

Usage:
    surface_flux.py --tube-r=<r>
"""
import os
import sys
import glob

import numpy as np
import yt.mods as ytm

import pysac.io
import pysac.io.gdf_writer
import pysac.io.yt_fields
import pysac.analysis.tube3D.process_utils as util
import pysac.analysis.tube3D.tvtk_tube_functions as ttf

sys.path.append('../')
from scripts import sacconfig
cfg = sacconfig.SACConfig()

try:
    import docopt
except ImportError:
    from script.extern import docopt

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

arguments = docopt.docopt(__doc__, version='Surface Analysis 13/11/13')

def glob_files(tube_r, search):
    files = glob.glob(os.path.join(cfg.data_dir,tube_r,search))
    files.sort()
    return files

def path_join(filename):
    return os.path.join(os.path.join(cfg.data_dir,'%s/'%tube_r),filename)

#Read Wave Flux HDF5 files in using yt
timeseries = ytm.load(os.path.join(cfg.gdf_dir,"*{}_fwave_0*.gdf".format(cfg.str_exp_fac)))
ds = timeseries[0]

#==============================================================================
# Define some crap
#==============================================================================
top_cut = -5
cube_slice = np.s_[:,:,:top_cut]
x_slice = np.s_[:,:,:,:top_cut]
cg = ds.h.grids[0]
#nlines is the number of fieldlines used in the surface
n_lines = 100
#the line is the fieldline to use as "the line"
line_n = 25
#==============================================================================

tube_r = arguments['--tube-r']

#Read in the velocity vtk files
linesurf = glob_files(tube_r,'Fieldline_surface*')

#Read in the surface indicies for the line
save_index = np.load(path_join("LineVar_%s_%s_%s_%s__%s_index.npy"%(cfg.driver, cfg.str_period,
                                                          cfg.amp, tube_r, cfg.str_exp_fac)))
save_index = np.asarray(save_index, dtype=np.int)

#Divide up the time steps to each process
all_indices = np.arange(0, save_index.shape[0], 1, dtype=int)
if mpi:
    chunks = np.array_split(all_indices, size)
    rank_indices = comm.scatter(chunks, root=0)
else:
    rank_indices = all_indices
print rank_indices

#Create line output arrays
#Number of indicies for this process:
number_ind = len(rank_indices)

Fpar_line = np.zeros([number_ind, save_index.shape[1]])
Fperp_line = np.zeros([number_ind, save_index.shape[1]])
Fphi_line = np.zeros([number_ind, save_index.shape[1]])

#Iterate over the timeseries to save the waveflux surfaces
#    for i, ds in enumerate(timeseries):
#for k,i in enumerate(rank_indices):
def do_step(k,i):
    print i
    ds = timeseries[i]

    #Get the surface and scalars from the first one in the series
    surf_poly = ttf.read_step(linesurf[i])

    normals = ttf.get_data(surf_poly, 'perp')
    parallels = ttf.get_data(surf_poly, 'par')
    torsionals = ttf.get_data(surf_poly, 'phi')

    #Extract the wave flux mayavi field from the yt dataset
    fwfield = util.get_mlab_field_yt(ds, 'wave_flux_x', 'wave_flux_y', 'wave_flux_z', cube_slice=cube_slice)

    #Interpolate the Wave flux to the velocity surface
    surface_fwave_filter, surface_fwave = ttf.interpolate_vectors(fwfield.outputs[0],
                                                                 surf_poly)

    Fwperp, Fwpar, Fwphi = ttf.get_surface_velocity_comp(surface_fwave,
                                                         normals, torsionals,
                                                         parallels)

    ttf.write_wave_flux(path_join("WaveFlux_%s_%s_%s_%05i.vtp"%(cfg.driver, cfg.str_period, cfg.amp, i+1)),
                             surf_poly, parallels, normals, torsionals, Fwpar, Fwperp, Fwphi)

    Fpar_line[k] = Fwpar[save_index[i]]
    Fperp_line[k] = Fwperp[save_index[i]]
    Fphi_line[k] = Fwphi[save_index[i]]

#Put the loop in as a function call to fix a memory leak!
for k,i in enumerate(rank_indices):
    do_step(k,i)

if mpi:
    comm.Barrier()
    Fpar_line_r0 = comm.gather(Fpar_line, root=0)
    Fperp_line_r0 = comm.gather(Fperp_line, root=0)
    Fphi_line_r0 = comm.gather(Fphi_line, root=0)
else:
    #Make a bonus leading axis so it looks like the MPI arrays
    Fpar_line_r0 = Fpar_line[None]
    Fperp_line_r0 = Fperp_line[None]
    Fphi_line_r0 = Fphi_line[None]
if rank == 0:
    Fpar_line = np.concatenate(Fpar_line_r0,axis=0)
    Fperp_line = np.concatenate(Fperp_line_r0,axis=0)
    Fphi_line = np.concatenate(Fphi_line_r0,axis=0)

    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_Fpar.npy"%(cfg.driver, cfg.str_period,
                                                              cfg.amp, tube_r, cfg.str_exp_fac)),Fpar_line)
    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_Fperp.npy"%(cfg.driver, cfg.str_period,
                                                              cfg.amp, tube_r, cfg.str_exp_fac)),Fperp_line)
    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_Fphi.npy"%(cfg.driver, cfg.str_period,
                                                              cfg.amp, tube_r, cfg.str_exp_fac)),Fphi_line)
