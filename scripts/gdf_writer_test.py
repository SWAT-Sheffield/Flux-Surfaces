# -*- coding: utf-8 -*-

import h5py
import numpy as np
import astropy.units as u

from pysac.io.legacy import VACfile
from pysac.io.legacy.gdf_converter import convert_w_3D, write_gdf
from pysac.io.legacy.util import mag_convert

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.rank


input_fname = '/fastdata/smq11sjm/inidata/3D_tube_128_128_128.ini'

vfile = VACfile(input_fname)

mu = 1.25663706e-6
#==============================================================================
# Monkey patch varnames and make fields dictionary
#==============================================================================
def get_header_fields(vfile):
    header = vfile.header
#    header['varnames'] = cfg.varnames
#    indices = range(0,len(header['varnames']))
#    w_ = dict(zip(header['varnames'],indices))
#    
    #Convert w
    w = mag_convert(vfile.w, vfile.w_)
    fields = convert_w_3D(np.ascontiguousarray(w), vfile.w_)
    
    for field in fields.values():
        unit = field['field'].unit
        x = np.ascontiguousarray(np.rollaxis(np.array(field['field']),0,3))
        field['field'] = u.Quantity(x, unit=unit)
    return header, fields

header, fields = get_header_fields(vfile)

#==============================================================================
# Decide how big the whole array is and what slice of it this processor has
#==============================================================================
nx = header['nx']
#Split sizes
n0 = 2
n1 = 2
n2 = 4

full_nx = [nx[0]*n0, nx[1]*n1, nx[2]*n2]
header['nx'] = full_nx
print rank, nx, full_nx, vfile.x.shape
coords = np.zeros(3)
if rank < n0:
    coords[0] = rank
elif (rank == n0):
    coords[1] = 1 
elif rank < (n0 + n1):
    coords[0] = rank - n0
    coords[1] = rank - n1
else:
    coords[2] = rank / (n0 * n1)
    rank2 = rank - (coords[2] * n2)
    
    if rank2 < n0:
        coords[0] = rank2
    elif rank2 == n0:
        coords[1] = 1
    elif rank2 < (n0 + n0):
        coords[0] = rank2 - n0
        coords[1] = rank2 - n1

s = map(int, coords * nx)
e = map(int, (coords+1) * nx)
arr_slice =  np.s_[s[1]:e[1],s[2]:e[2],s[0]:e[0]]

#==============================================================================
# Reconstruct the whole x array on all the processors
#==============================================================================
x_slice =  np.s_[:,s[0]:e[0],s[1]:e[1],s[2]:e[2]]
x = np.zeros([3]+full_nx)
x[x_slice] = vfile.x

x_g = comm.gather(x, root=0)
if rank == 0:
    x_0 = np.sum(x_g,axis=0)
    x_xyz = np.zeros(x_0.shape)
    x_xyz[0] = np.rollaxis(x_0[1],0,3)#x
    x_xyz[1] = np.rollaxis(x_0[2],0,3)#y
    x_xyz[2] = np.rollaxis(x_0[0],0,3)#z
else:
    x_xyz = None

x = comm.bcast(x_xyz, root=0)
x *= u.meter

#==============================================================================
# Save a gdf file
#==============================================================================
for i in range(0,vfile.num_records,1):
    print '#', rank, '#', "read step %i"%i
    vfile.read_timestep(i+1)
    header, fields = get_header_fields(vfile)
    header['nx'] = full_nx

    f = h5py.File('/archive/gdf_testing/test'+'_%05i.gdf'%(i+1),
                  'w', driver='mpio', comm=MPI.COMM_WORLD)

    write_gdf(f, header, x, fields, arr_slice=arr_slice,
              data_author = "Stuart Mumford",
              data_comment = "Converted from outfiles in parallel")
print "rank %i finishes"%rank
