#!/usr/bin/env python2
"""
Configure repository, configure and compile SAC and analysis

Usage:
    configure.py set SAC [--compiler=<compiler>] [--compiler_flags=<flg>] [--vac_modules=<mod>] [--runtime=<s>] [--mpi_config=<cfg>] [--varnames=<s>]
    configure.py set driver [--driver=<driver>] [--period=<s>] [--exp_fac=<exp>] [--amp=<amp_str>] [--fort_amp=<fort_amp>]
    configure.py set analysis [--tube_radii=<radii>]
    configure.py set data [--ini_dir=<dir>] [--out_dir=<dir>] [--data_dir=<dir>] [--gdf_dir=<dir>]
    configure.py print [<section>]
    configure.py compile SAC [--clean]

Options:
    --compiler=COMPILER  The FORTRAN compiler to use to compile SAC
    --compiler_flags=FLAGS  Arguments to pass to the compiler
    --vac_modules=MODS  VAC modules to include in compilation
    --runtime=S  Physical time runtime of the simulation (in seconds)
    --mpi_config=CFG  The mpi conifg to use
    --varnames=NAMES  A list of physical varnames for all the varibles in w
    --driver=DRI  The name of the driver
    --period=S  The driver period in seconds
    --exp_fac=E  The expqansion factor of the driver
    --amplitude=AMP  String representation of the driver ampliude
    --fort_amp=FAMP  FORTRAN code to calculate the driver amplitude
    --tube-radii=R  A list of flux surface radii to compute analysis for
    --ini_dir=DIR The data directory where the initial conditions are located.
    --out_dir=DIR  SAC output data directory
    --data_dir=DIR  ?
    --gdf_dir=DIR  The dir for GDF output, identifier will be appended as a dir
    --clean  Clean SAC compile directory before compiling
"""
import os
import sys

try:
    import docopt
except ImportError:
    import scripts.extern.docopt as docopt
arguments = docopt.docopt(__doc__, version='SAC Configuration')

from scripts import sacconfig

#==============================================================================
# User Configuration
#==============================================================================
cfg = sacconfig.SACConfig()

if arguments['set']:
    #Strip out all optional arguments the user has set
    opts = {}
    for k,v in arguments.items():
        if k[:2] == '--' and v is not None:
            opts.update({k[2:]:v})

    #Update the config with the options
    for k,v in opts.items():
        cfg.__setattr__(k, v)

    #Save the config file
    cfg.save_cfg()

if arguments['print']:
    if not arguments['<section>']:
        section = 'all'
    else:
        section = arguments['<section>']

    cfg.print_config(section)

driver = cfg.driver
period = cfg.period
post_amp = cfg.amp
exp_fac = cfg.str_exp_fac
out_dir = cfg.out_dir
tube_radii = cfg.tube_radii

identifier = cfg.get_identifier()
job_name = "sac_%s_%s"%(driver,period)

def sac_path(path):
    return os.path.join(os.path.realpath('./sac/sac/'), path)

def check_file(path):
    # Make sure the current vac.par is a symlink
    path = sac_path(path)
    if os.path.isfile(path) or os.path.islink(path):
        if os.path.islink(path):
            os.remove(path)
        else:
            raise TypeError("{} is not a symlink, cannot safely remove".format(path))

if arguments['compile'] and arguments['SAC']:
    # Get the template vac.par into place
    check_file('vac.par')
    os.symlink(os.path.realpath('./scripts/vac_config.par'),
               os.path.realpath(sac_path('vac.par')))

    check_file('src/vacusrpar.t.Slog')
    os.symlink(os.path.realpath('./scripts/vacusrpar.t.Slog'),
               sac_path('src/vacusrpar.t.Slog'))

    check_file('src/vacusr.t.Slog')
    os.symlink(os.path.realpath('./scripts/vacusr.t.Slog'),
               sac_path('src/vacusr.t.Slog'))

    #==============================================================================
    # Process vac.par
    #==============================================================================
    f = open(sac_path('vac.par'), 'r')
    f_lines = f.readlines()

    for i,line in enumerate(f_lines):
        if line.strip().startswith("filenameini="):
            f_lines[i] = '\t' + "filenameini='" + os.path.join(cfg.ini_dir,
                                        "3D_tube_128_128_128.ini'\n")
        if line.strip().startswith("filename="):
            f_lines[i] = '\t' + "filename='" + os.path.join(out_dir,
                                    "3D_tube128_%s.log'\n"%identifier)
            f_lines[i+1] = "\t\t'" + os.path.join(out_dir,
                                    "3D_tube128_%s.out'\n"%identifier)
        if line.strip().startswith("wnames="):
            f_lines[i] = '\t' + "wnames='" + ' '.join(cfg.varnames).encode('ascii') +"'\n"
    f.close()

    for i in range(len(f_lines)):
        f_lines[i] = f_lines[i].encode('ascii')

    #Truncate and overwrite
    f = open(sac_path('vac.par'), 'w')
    f.writelines(f_lines)
    f.close()

    #==============================================================================
    # Process vacusr.t.Slog
    #==============================================================================
    f = open(sac_path('src/vacusr.t.Slog'), 'r')
    f_lines = f.readlines()

    for i,line in enumerate(f_lines):
        if line.strip() == "!### AMPLITUDE ###":
            f_lines[i+1] = '    AA = %s\n'%cfg.fort_amp
        if line.strip() == "!### PERIOD ###":
            f_lines[i+1] = '    s_period = %s\n'%cfg.period
        if line.strip() == "!### EXP_FAC ###":
            f_lines[i+1] = '    B = %s\n'%cfg.exp_fac
    f.close()

    #Truncate and overwrite
    f = open(sac_path('src/vacusr.t.Slog'), 'w')
    f.writelines(f_lines)
    f.close()

    #==============================================================================
    # Compile Things
    #==============================================================================
    #Compile VAC
    os.chdir(sac_path("src"))
    os.system('./setvac -u=Slog -p=mhd')
    if arguments['--clean']:
        os.system('./sac_fabricate.py --clean')
    else:
        os.system('./sac_fabricate.py')
    os.chdir(os.path.realpath('.'))
