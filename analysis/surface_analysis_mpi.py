#!/usr/bin/env python2
"""
Create a flux surface and decompose varibales onto it, in MPI parallel or serial

Usage:
    surface_analysis_mpi.py --tube-r=<r>
"""
import os
import sys
import glob

import h5py
import numpy as np
import astropy.units as u
from tvtk.api import tvtk
from mayavi.tools.sources import vector_field
import yt

import pysac.io
import pysac.yt
import pysac.io.gdf_writer
import pysac.analysis.tube3D.process_utils as util
import pysac.analysis.tube3D.tvtk_tube_functions as ttf

sys.path.append('../')
from scripts import sacconfig

try:
    import docopt
except ImportError:
    from scripts.extern import docopt

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
arguments = docopt.docopt(__doc__, version='Surface Analysis 13/11/13')
tube_r = arguments['--tube-r']

cfg = sacconfig.SACConfig()

#Make the pull paths
data_dir = os.path.join(cfg.data_dir,'%s/'%tube_r)
gdf_files = glob.glob(os.path.join(cfg.gdf_dir, cfg.get_identifier()+'_0*.gdf'))
gdf_files.sort()

if rank == 0:
    print "Configuration:"
    print 'Identifier:', cfg.get_identifier()
    print 'tube_r:', tube_r
    print 'data_dir:', data_dir
    print 'gdf_dir:', cfg.gdf_dir

if rank == 0: #Prevents race condition where one processes creates the dir
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

def path_join(filename):
    return os.path.join(data_dir,filename)

#==============================================================================
# TimeSeries, data slices and other constants
#==============================================================================
#Create the yt time series
#ts = yt.load(gdf_files)
np.save(os.path.join(cfg.data_dir,'Times_{}.npy'.format(cfg.get_identifier())),
        [pysac.yt.SACGDFDataset(fname).current_time for fname in gdf_files])
#Define a var to limit iterations, no limt = len(ts)
max_n = 16#len(ts)

top_cut = -5
cube_slice = np.s_[:,:,:top_cut]
x_slice = np.s_[:,:,:,:top_cut]
ds = pysac.yt.SACGDFDataset(gdf_files[0])
cg = ds.index.grids[0]
xmax, ymax, zmax = cg['density_pert'][cube_slice].shape #- 1# ts[0].domain_dimensions - 1
#nlines is the number of fieldlines used in the surface
n_lines = 100
#the line is the fieldline to use as "the line"
line_n = 25

#==============================================================================
# Do Serial compute of seed points
# Use h5py as it is very slightly faster than yt(2.6):
# https://hub.yt-project.org/nb/jlrbq9
#==============================================================================
if rank == 0:
    f = h5py.File(gdf_files[0])

    save_times = []
    xmax, ymax, zmax = ds.domain_dimensions - 1
    zmax += top_cut
    #nlines is the number of fieldlines used in the surface
    n_lines = 100
    #the line is the fieldline to use as "the line"
    line_n = 25

    seeds_slice = np.s_[:,:,top_cut-1]#slicing upto -5 is not the same as indexing -5

    vfield = vector_field(f['data/grid_0000000000']['velocity_x'][seeds_slice] / 1e3 * 1e-2,
                          f['data/grid_0000000000']['velocity_y'][seeds_slice] / 1e3 * 1e-2,
                          f['data/grid_0000000000']['velocity_z'][seeds_slice] / 1e3 * 1e-2,
                          name="Velocity Field", figure=None).outputs[0]

    #Define domain parameters
    domain = {'xmax':xmax,'ymax':ymax,'zmax':0}
    #Create initial seed points in tvtk.PolyData
    surf_seeds_poly = ttf.make_circle_seeds(n_lines, int(tube_r[1:]), **domain)
    #Make surface using seeds, surf filter and contour
    #Extract all times
    next_seeds = np.array(surf_seeds_poly.points)
    next_seeds[:,-1] = zmax
    surf_seeds = [next_seeds]
    t1 = f['simulation_parameters'].attrs['current_time']
    save_times.append(t1)

    for i in range(0,max_n):
        print i
        f = h5py.File(gdf_files[i])
        t1 = save_times[-1]
        t2 = f['simulation_parameters'].attrs['current_time']
        save_times.append(t2)
        vfield = vector_field(f['data/grid_0000000000']['velocity_x'][seeds_slice] / 1e3 * 1e-2,
                              f['data/grid_0000000000']['velocity_y'][seeds_slice] / 1e3 * 1e-2,
                              f['data/grid_0000000000']['velocity_z'][seeds_slice] / 1e3 * 1e-2,
                              name="Velocity Field", figure=None).outputs[0]

        next_seeds = ttf.move_seeds(surf_seeds_poly, vfield, t2-t1)
        next_seeds[:,-1] = zmax
        surf_seeds.append(next_seeds)
else:
    surf_seeds = None
    save_times = []

## tvtk objects are presumeably not pikleable, therefore cannot be trasmitted
## via MPI, create the tvtk object upon recieve.
surf_seeds_arr = np.array(surf_seeds)

if rank == 0:
    #Save the seeds for good measure
    np.save(path_join("SeedPoints_{}_{}.npy".format(cfg.get_identifier(), tube_r)), surf_seeds_arr)

save_times = np.array(save_times)
if mpi:
    #Give all processes the surf_seeds
    surf_seeds_arr = comm.bcast(surf_seeds_arr, root=0)
    save_times = comm.bcast(save_times, root=0)

surf_seeds = []
for seeds in surf_seeds_arr:
    pd = tvtk.PolyData()
    pd.points = seeds
    surf_seeds.append(pd)

#Divide up the time steps to each process
all_indices = np.arange(0,max_n,1,dtype=int)

if mpi:
    chunks = np.array_split(all_indices, size)
    rank_indices = comm.scatter(chunks, root=0)
else:
    rank_indices = all_indices

#==============================================================================
# Begin Parallel compute
#==============================================================================
ds = pysac.yt.SACGDFDataset(gdf_files[0])
#ds = ts[0]
#Read in all data from hdf5 file
bfield, vfield, density, valf, cs, beta = util.get_yt_mlab(ds, cube_slice, flux=True)

bpert = util.yt_to_mlab_vector(ds,
                               'mag_field_x_pert',
                               'mag_field_y_pert',
                               'mag_field_z_pert',
                               cube_slice,
                               field_name="Magnetic Pertubation")

f_wave = pysac.analysis.get_wave_flux_yt(ds)  # , B_to_SI=1e-4, V_to_SI=1e-2, Pk_to_SI=1e-1)

fwfield = vector_field(f_wave[0,:,:,:][cube_slice],
                       f_wave[1,:,:,:][cube_slice],
                       f_wave[2,:,:,:][cube_slice],
                       name="Wave Flux",figure=None)

#Make surface using seeds, surf filter and contour
surf_field_lines, surface = ttf.create_flux_surface(bfield.outputs[0], surf_seeds[0])

#Make the PolyDataNormals object
poly_norms = ttf.make_poly_norms(surface.output)

#Interpolate the vfield to the surface
surface_vel_filter, surface_velocities = ttf.interpolate_vectors(vfield.outputs[0], poly_norms.output)
surface_mag_filter, surface_bfield = ttf.interpolate_vectors(bfield.outputs[0], poly_norms.output)

surface_bpert_filter, surface_bpert = ttf.interpolate_vectors(bpert.outputs[0], poly_norms.output)

#Interpolate the Wave Flux
surface_fwave_filter, surface_fwave = ttf.interpolate_vectors(fwfield.outputs[0], poly_norms.output)

#Interpolate the vfield to the surface
surface_den_filter, surface_density = ttf.interpolate_vectors(density.outputs[0], poly_norms.output)
surface_va_filter, surface_va = ttf.interpolate_vectors(valf.outputs[0], poly_norms.output)
surface_cs_filter, surface_cs = ttf.interpolate_vectors(cs.outputs[0], poly_norms.output)
surface_beta_filter, surface_beta = ttf.interpolate_vectors(beta.outputs[0], poly_norms.output)

#Make the line
the_line = ttf.get_the_line(bfield.outputs[0], surf_seeds[0], line_n)

#Number of indicies for this process:
number_ind = len(rank_indices)

#==============================================================================
# Set up arrays to store all the varibles for all the timesteps
#==============================================================================
# Velocity and Surface Info
save_points = np.zeros([number_ind, len(the_line.output.points),3])
save_vpar = np.zeros([number_ind, len(the_line.output.points)])
save_vperp = np.zeros([number_ind, len(the_line.output.points)])
save_vphi = np.zeros([number_ind, len(the_line.output.points)])
save_index = np.zeros([number_ind, len(the_line.output.points)])

save_bpertpar = np.zeros([number_ind, len(the_line.output.points)])
save_bpertperp = np.zeros([number_ind, len(the_line.output.points)])
save_bpertphi = np.zeros([number_ind, len(the_line.output.points)])

#Flux Lines
Fpar_line = np.zeros([number_ind, len(the_line.output.points)])
Fperp_line = np.zeros([number_ind, len(the_line.output.points)])
Fphi_line = np.zeros([number_ind, len(the_line.output.points)])

# Scalar Vars
density_line = np.zeros([number_ind, len(the_line.output.points)])
va_line = np.zeros([number_ind, len(the_line.output.points)])
beta_line = np.zeros([number_ind, len(the_line.output.points)])
cs_line = np.zeros([number_ind, len(the_line.output.points)])

# Whole surface Avgs
Fpar_avg = np.zeros([number_ind])
Fperp_avg = np.zeros([number_ind])
Fphi_avg = np.zeros([number_ind])

for i,n in enumerate(rank_indices):
    ds = pysac.yt.SACGDFDataset(gdf_files[n])
    #ds = ts[n]
    [bfield, vfield, density,
     valf, cs, beta] = util.process_next_step_yt(ds, cube_slice, bfield, vfield,
                                                 density, valf, cs, beta)

    bpert = util.update_yt_to_mlab_vector(bpert, ds,
                                          'mag_field_x_pert',
                                          'mag_field_y_pert',
                                          'mag_field_z_pert',
                                          cube_slice)
    #Update surface
    ttf.update_flux_surface(surf_seeds[n], surf_field_lines, surface)
    #Get velocities at surface
    surface_velocities = ttf.update_interpolated_vectors(False,
                                                         surface_vel_filter)
    #Get surface magfield
    surface_bfield = ttf.update_interpolated_vectors(False, surface_mag_filter)

    surface_bpert = ttf.update_interpolated_vectors(False, surface_bpert_filter)

    #Get vectors at surface
    normals, torsionals, parallels = ttf.get_surface_vectors(poly_norms,
                                                         surface_bfield)
    #Get velocity componets at surface
    vperp, vpar, vphi = ttf.get_surface_velocity_comp(surface_velocities,
                                                      normals, torsionals,
                                                      parallels)

    #Get magnetic field pertubation componets at surface
    bpertperp, bpertpar, bpertphi = ttf.get_surface_velocity_comp(surface_bpert,
                                                           normals, torsionals,
                                                           parallels)

    surface_density = ttf.update_interpolated_scalars(surface.output, surface_den_filter)
    surface_va = ttf.update_interpolated_scalars(surface.output, surface_va_filter)
    surface_beta = ttf.update_interpolated_scalars(surface.output, surface_beta_filter)
    surface_cs = ttf.update_interpolated_scalars(surface.output, surface_cs_filter)

    #Save to file
    writer = ttf.PolyDataWriter(path_join("Fieldline_surface_{}_{:05d}.vtp".format(cfg.get_identifier(), n+1)),
                                surface.output)
    writer.add_array(perp=normals,
                     par=parallels,
                     phi=torsionals,
                     vperp=vperp,
                     vpar=vpar,
                     bpertperp=bpertperp,
                     bpertpar=bpertpar,
                     bpertphi=bpertphi,
                     surface_density=surface_density,
                     surface_va=surface_va,
                     surface_beta=surface_beta,
                     surface_cs=surface_cs)
    writer.write()

#==============================================================================
# Save out a long one field line
#==============================================================================
    the_line = ttf.update_the_line(the_line, surf_seeds[n].points, line_n,
                                   len(the_line.output.points))

    #line varibles
    surf_line_index, surf_line_points = ttf.get_surface_indexes(surface.output,
                                                                the_line)

    save_index[i] = surf_line_index
    save_points[i] = surf_line_points

    save_vpar[i] = vpar[surf_line_index]
    save_vperp[i] = vperp[surf_line_index]
    save_vphi[i] = vphi[surf_line_index]

    save_bpertpar[i] = bpertpar[surf_line_index]
    save_bpertperp[i] = bpertperp[surf_line_index]
    save_bpertphi[i] = bpertphi[surf_line_index]

    # Scalar Fluxes
    density_line[i] = surface_density[surf_line_index]
    va_line[i] = surface_va[surf_line_index]
    beta_line[i] = surface_beta[surf_line_index]
    cs_line[i] = surface_cs[surf_line_index]

#==============================================================================
# Wave Flux
#==============================================================================
    f_wave = pysac.analysis.get_wave_flux_yt(ds)  # , B_to_SI=1e-4, V_to_SI=1e-2, Pk_to_SI=1e-1)

    fwfield.set(vector_data = np.rollaxis(np.array([f_wave[0,:,:,:][cube_slice],
                                                    f_wave[1,:,:,:][cube_slice],
                                                    f_wave[2,:,:,:][cube_slice]]),
                                          0, 4))

    surface_fwave = ttf.update_interpolated_vectors(False, surface_fwave_filter)

    Fwperp, Fwpar, Fwphi = ttf.get_surface_velocity_comp(surface_fwave,
                                                         normals, torsionals,
                                                         parallels)

    # Save out the Wave Flux to a smaller file, with vector info
    writer = ttf.PolyDataWriter(path_join("WaveFlux_{}_{:05d}.vtp".format(cfg.get_identifier(), n+1)),
                                surface.output)
    writer.add_array(par=parallels,
                     perp=normals,
                     phi=torsionals,
                     Fwpar=Fwpar,
                     Fwperp=Fwperp,
                     Fwphi=Fwphi)
    writer.write()

    # Line Flux
    Fpar_line[i] = Fwpar[surf_line_index]
    Fperp_line[i] = Fwperp[surf_line_index]
    Fphi_line[i] = Fwphi[surf_line_index]

    # Whole Surface Avg
    Fperp_avg[i] = np.mean(Fwperp)
    Fpar_avg[i] = np.mean(Fwpar)
    Fphi_avg[i] = np.mean(Fwphi)

    #Save waveflux to the gdf file:
    fname, ext = os.path.splitext(ds.parameter_filename)
    fprefix = fname[:-6]
    fnumber = fname[-6:]
    new_name = fprefix + '_fwave' + fnumber + ext

    f = h5py.File(new_name, mode='w')
    simulation_params = pysac.io.gdf_writer.SimulationParameters()
    simulation_params['dimensionality'] = 3
    simulation_params['domain_dimensions'] = ds.domain_dimensions
    simulation_params['current_time'] = ds.current_time
    simulation_params['domain_left_edge'] = ds.domain_left_edge
    simulation_params['domain_right_edge'] = ds.domain_right_edge
    simulation_params['num_ghost_zones'] = [0]
    simulation_params['field_ordering'] = 0
    simulation_params['boundary_conditions'] = np.zeros([6], dtype=int)+2

    f = pysac.io.gdf_writer.create_file(f, simulation_params, ds.domain_dimensions,
                                        data_author="Flux-Surfaces")

    fwave_x = u.Quantity(f_wave[0,:,:,:], unit='W/m2')
    pysac.io.gdf_writer.write_field(f, fwave_x, 'wave_flux_x', 'x Compoment of Wave Energy Flux')
    fwave_y = u.Quantity(f_wave[1,:,:,:], unit='W/m2')
    pysac.io.gdf_writer.write_field(f, fwave_y, 'wave_flux_y', 'y Compoment of Wave Energy Flux')
    fwave_z = u.Quantity(f_wave[2,:,:,:], unit='W/m2')
    pysac.io.gdf_writer.write_field(f, fwave_z, 'wave_flux_z', 'z Compoment of Wave Energy Flux')
    f.close()

    print "%i: Done Flux %i \n%i: Step %i / %i"%(rank,n,rank,i,len(rank_indices))

#==============================================================================
# Gather up the single line arrays and save to a file
#==============================================================================
print "%i: Finishing up"%rank

if mpi:
    comm.Barrier()
    #Gather the data
    save_points_r0 = comm.gather(save_points, root=0)
    save_index_r0= comm.gather(save_index, root=0)

    save_vpar_r0 = comm.gather(save_vpar, root=0)
    save_vperp_r0 = comm.gather(save_vperp, root=0)
    save_vphi_r0 = comm.gather(save_vphi, root=0)

    save_bpertpar_r0 = comm.gather(save_bpertpar, root=0)
    save_bpertperp_r0 = comm.gather(save_bpertperp, root=0)
    save_bpertphi_r0 = comm.gather(save_bpertphi, root=0)

    Fpar_line_r0 = comm.gather(Fpar_line, root=0)
    Fperp_line_r0 = comm.gather(Fperp_line, root=0)
    Fphi_line_r0 = comm.gather(Fphi_line, root=0)

    density_line_r0 = comm.gather(density_line, root=0)
    va_line_r0 = comm.gather(va_line, root=0)
    beta_line_r0 = comm.gather(beta_line, root=0)
    cs_line_r0 = comm.gather(cs_line, root=0)

    Fperp_avg_r0 = comm.gather(Fwperp, root=0)
    Fpar_avg_r0 = comm.gather(Fwpar, root=0)
    Fphi_avg_r0 = comm.gather(Fwphi, root=0)
else:
    #Make a bonus leading axis so it looks like the MPI arrays
    save_points_r0 = save_points[None]
    save_index_r0= save_index[None]

    save_vpar_r0 = save_vpar[None]
    save_vperp_r0 = save_vperp[None]
    save_vphi_r0 = save_vphi[None]

    save_bpertpar_r0 = save_bpertpar[None]
    save_bpertperp_r0 = save_bpertperp[None]
    save_bpertphi_r0 = save_bpertphi[None]

    Fpar_line_r0 = Fpar_line[None]
    Fperp_line_r0 = Fperp_line[None]
    Fphi_line_r0 = Fphi_line[None]

    Fpar_avg_r0 = Fpar_avg[None]
    Fperp_avg_r0 = Fperp_avg[None]
    Fphi_avg_r0 = Fphi_avg[None]

    density_line_r0 = density_line[None]
    va_line_r0 = va_line[None]
    beta_line_r0 = beta_line[None]
    cs_line_r0 = cs_line[None]

if rank == 0:
    save_points = np.concatenate(save_points_r0, axis=0)
    save_index = np.concatenate(save_index_r0, axis=0)

    save_vpar = np.concatenate(save_vpar_r0, axis=0)
    save_vperp = np.concatenate(save_vperp_r0, axis=0)
    save_vphi = np.concatenate(save_vphi_r0, axis=0)

    save_bpertpar = np.concatenate(save_bpertpar_r0, axis=0)
    save_bpertperp = np.concatenate(save_bpertperp_r0, axis=0)
    save_bpertphi = np.concatenate(save_bpertphi_r0, axis=0)

    Fpar_line = np.concatenate(Fpar_line_r0, axis=0)
    Fperp_line = np.concatenate(Fperp_line_r0, axis=0)
    Fphi_line = np.concatenate(Fphi_line_r0, axis=0)

    Fpar_avg = np.concatenate(Fpar_avg_r0, axis=0)
    Fperp_avg = np.concatenate(Fperp_avg_r0, axis=0)
    Fphi_avg = np.concatenate(Fphi_avg_r0, axis=0)

    density_line = np.concatenate(density_line_r0, axis=0)
    va_line = np.concatenate(va_line_r0, axis=0)
    beta_line = np.concatenate(beta_line_r0, axis=0)
    cs_line = np.concatenate(cs_line_r0, axis=0)

    np.save(path_join("LineVar_{}_points.npy".format(cfg.get_identifier())), save_points)
    np.save(path_join("LineVar_{}_times.npy".format(cfg.get_identifier())), save_times)
    np.save(path_join("LineVar_{}_index.npy".format(cfg.get_identifier())), save_index)
    
    np.save(path_join("LineVar_{}_vphi.npy".format(cfg.get_identifier())), save_vphi)
    np.save(path_join("LineVar_{}_vperp.npy".format(cfg.get_identifier())), save_vperp)
    np.save(path_join("LineVar_{}_vpar.npy".format(cfg.get_identifier())), save_vpar)
    
    np.save(path_join("LineVar_{}_bpertphi.npy".format(cfg.get_identifier())), save_bpertphi)
    np.save(path_join("LineVar_{}_bpertperp.npy".format(cfg.get_identifier())), save_bpertperp)
    np.save(path_join("LineVar_{}_bpertpar.npy".format(cfg.get_identifier())), save_bpertpar)
    
    np.save(path_join("LineVar_{}_rho.npy".format(cfg.get_identifier())), density_line)
    np.save(path_join("LineVar_{}_va.npy".format(cfg.get_identifier())), va_line)
    np.save(path_join("LineVar_{}_cs.npy".format(cfg.get_identifier())), cs_line)
    np.save(path_join("LineVar_{}_beta.npy".format(cfg.get_identifier())), beta_line)
    
    np.save(path_join("LineFlux_{}_Fpar.npy".format(cfg.get_identifier())), Fpar_line)
    np.save(path_join("LineFlux_{}_Fperp.npy".format(cfg.get_identifier())), Fperp_line)
    np.save(path_join("LineFlux_{}_Fphi.npy".format(cfg.get_identifier())), Fphi_line)
    
    np.save(path_join("AverageFlux_{}_Fperp.npy".format(cfg.get_identifier())), Fperp_avg)
    np.save(path_join("AverageFlux_{}_Fpar.npy".format(cfg.get_identifier())), Fpar_avg)
    np.save(path_join("AverageFlux_{}_Fphi.npy".format(cfg.get_identifier())), Fphi_avg)
