# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 16:45:11 2014

Calculate and save to disk the surface Flux averages

Usage:
    whole_surface_avgs.py
"""
import os
import sys
import glob

import numpy as np

import pysac.analysis.tube3D.tvtk_tube_functions as ttf

sys.path.append('../')
from scripts import sacconfig
cfg = sacconfig.SACConfig()

try:
    import docopt
except ImportError:
    from script.extern import docopt

arguments = docopt.docopt(__doc__, version='Surface Analysis 13/11/13')

driver = cfg.driver
post_amp = cfg.amp
period = cfg.str_period

#Add the '_' to exp_fac
if cfg.exp_fac:
    exp_fac = '_' + cfg.str_exp_fac
else:
    exp_fac=''

def glob_files(tube_r, search):
    files = glob.glob(os.path.join(cfg.data_dir,tube_r,search))
    files.sort()
    return files

def path_join(filename):
    return os.path.join(os.path.join(cfg.data_dir,'%s/'%tube_r),filename)

for tube_r in cfg.tube_radii:
    print tube_r
    wave = glob_files(tube_r, 'Wave*')
    
    Fpar = np.zeros([len(wave)])
    Fperp = np.zeros([len(wave)])
    Fphi = np.zeros([len(wave)])
    
    for i, awave in enumerate(wave):
        surf_poly = ttf.read_step(awave)
    
        Fwperp = ttf.get_data(surf_poly, 'Fwperp')
        Fwpar = ttf.get_data(surf_poly, 'Fwpar')
        Fwphi = ttf.get_data(surf_poly, 'Fwphi')
        
        Fperp[i] = np.mean(Fwperp)
        Fpar[i] = np.mean(Fwpar)
        Fphi[i] = np.mean(Fwphi)
        
        Fwpar[np.abs(Fwpar)<1e-5], Fwperp[np.abs(Fwperp)<1e-5], Fwphi[np.abs(Fwphi)<1e-5] = 0., 0., 0.

        Fwtot = np.sqrt(Fwpar**2 + Fwperp**2 + Fwphi**2)
        Fwpar, Fwperp, Fwphi = (Fwpar/Fwtot)*100, (Fwperp/Fwtot)*100, (Fwphi/Fwtot)*100
        
        Fwpar, Fwperp, Fwphi = map(np.abs, [Fwpar, Fwperp, Fwphi])
        
        Fwpar = np.mean(np.ma.masked_array(Fwpar,np.isnan(Fwpar)))
        Fwperp = np.mean(np.ma.masked_array(Fwperp,np.isnan(Fwperp)))
        Fwphi = np.mean(np.ma.masked_array(Fwphi,np.isnan(Fwphi)))
    
    np.save(path_join("AveragePercentFlux_%s_%s_%s_%s_%s_Fperp.npy"%(cfg.driver, cfg.str_period,
                                                                  cfg.amp, tube_r, cfg.str_exp_fac)),Fwperp)
    np.save(path_join("AveragePercentFlux_%s_%s_%s_%s_%s_Fpar.npy"%(cfg.driver, cfg.str_period,
                                                                  cfg.amp, tube_r, cfg.str_exp_fac)),Fwpar)
    np.save(path_join("AveragePercentFlux_%s_%s_%s_%s_%s_Fphi.npy"%(cfg.driver, cfg.str_period,
                                                                  cfg.amp, tube_r, cfg.str_exp_fac)),Fwphi)