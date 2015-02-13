# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 15:12:13 2014

@author: stuart
"""

import os
import glob

import numpy as np
import yt.mods as ytm

from sacconfig import SACConfig
cfg = SACConfig()


gdf_files = glob.glob(os.path.join(cfg.gdf_dir, cfg.get_identifier()+'_0*.gdf'))
gdf_files.sort()
ts = ytm.load(gdf_files)

np.save(os.path.join(cfg.data_dir,'Times_{}.npy'.format(cfg.get_identifier())), [ds.current_time for ds in ts])