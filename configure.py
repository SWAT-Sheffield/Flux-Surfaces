#!/usr/bin/env python2
"""
Configure repository, configure and compile SAC and analysis

Usage:
    configure.py set SAC [--compiler=<compiler>] [--compiler_flags=<flg>] [--vac_modules=<mod>] [--runtime=<s>] [--mpi_config=<cfg>] [--varnames=<s>] [--grid_size=<grid_size>] [--usr_script=<usr>]
    configure.py set driver [--period=<s>] [--exp_fac=<exp>] [--amp=<amp_str>] [--fort_amp=<fort_amp>] [--delta_x=<delta_x>] [--delta_y=<delta_y>] [--delta_z=<delta_z>] [--identifier=<identifier>]
    configure.py set analysis [--tube_radii=<radii>]
    configure.py set data [--ini_dir=<dir>] [--out_dir=<dir>] [--data_dir=<dir>] [--gdf_dir=<dir>]
    configure.py get <key>
    configure.py print [<section>]
    configure.py compile SAC [--clean]

Options:
    --compiler=COMPILER  The FORTRAN compiler to use to compile SAC
    --compiler_flags=FLAGS  Arguments to pass to the compiler
    --vac_modules=MODS  VAC modules to include in compilation
    --runtime=S  Physical time runtime of the simulation (in seconds)
    --mpi_config=CFG  The mpi conifg to use
    --varnames=NAMES  A list of physical varnames for all the varibles in w
    --grid_size=SIZE The grid size for each MPI processor including boundaries if fullgridini is True
    --usr_script=USR The vacusr.t suffix, where the scripts are in scripts/
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
import glob

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

if arguments['get']:
    print getattr(cfg, arguments['<key>'])

if arguments['print']:
    if not arguments['<section>']:
        section = 'all'
    else:
        section = arguments['<section>']

    cfg.print_config(section)

driver = cfg.usr_script
period = cfg.period
post_amp = cfg.amp
exp_fac = cfg.str_exp_fac
out_dir = cfg.out_dir
tube_radii = cfg.tube_radii

identifier = cfg.get_identifier()
job_name = "sac_%s_%s"%(driver,period)

# Root path of repo
root_path = os.path.realpath(os.path.dirname(__file__))

def sac_path(path):
    return os.path.join(root_path, 'sac/sac/', path)

def check_file(path):
    # Make sure the current vac.par is a symlink
    path = sac_path(path)
    if os.path.isfile(path) or os.path.islink(path):
        if os.path.islink(path):
            os.remove(path)
        else:
            raise TypeError("{0} is not a symlink, cannot safely remove".format(path))

if arguments['compile'] and arguments['SAC']:
    # Get the template vac.par into place
    check_file('vac.par')
    os.symlink(os.path.realpath('./scripts/vac_config.par'),
               os.path.realpath(sac_path('vac.par')))

    check_file('src/vacusrpar.t.{0}'.format(cfg.usr_script))
    os.symlink(os.path.realpath('./scripts/vacusrpar.t.{0}'.format(cfg.usr_script)),
               sac_path('src/vacusrpar.t.{0}'.format(cfg.usr_script)))

    check_file('src/vacusr.t.{0}'.format(cfg.usr_script))
    os.symlink(os.path.realpath('./scripts/vacusr.t.{0}'.format(cfg.usr_script)),
               sac_path('src/vacusr.t.{0}'.format(cfg.usr_script)))

    # Distribute ini file:
    os.chdir(cfg.ini_dir)
    if not os.path.isfile('3D_tube_128_128_128.ini'):
        raise ValueError("No initial conditions found, please download with ./run.py download ini")
    else:
        if not len(glob.glob('3D_tube_128_128_128_{0}_*.ini'.format(cfg.mpi_config))) == 16:
            print "Distributing ini file..."
            os.system('{0} 3D_tube_128_128_128.ini 3D_tube_128_128_128_{0}.ini'.format(sac_path('distribution'),
                                                                                 cfg.mpi_config))
    os.chdir(root_path)
    #==============================================================================
    # Process vac.par
    #==============================================================================
    f = open(sac_path('vac.par'), 'r')
    f_lines = f.readlines()

    for i,line in enumerate(f_lines):
        if line.strip().startswith("filenameini="):
            f_lines[i] = '    ' + "filenameini='" + os.path.join(cfg.ini_dir,
                                        "3D_tube_128_128_128_{0}.ini'\n".format(cfg.mpi_config))
        if line.strip().startswith("filename="):
            f_lines[i] = '    ' + "filename='" + os.path.join(out_dir,
                                    "3D_tube128_%s.log'\n"%identifier)
            f_lines[i+1] = "             '" + os.path.join(out_dir,
                                    "3D_tube128_%s.out'\n"%identifier)
        if line.strip().startswith("wnames="):
            f_lines[i] = '    ' + "wnames='" + ' '.join(cfg.varnames).encode('ascii') +"'\n"
        if line.strip().startswith("tmax="):
            f_lines[i] = '    ' + "tmax={0}\n".format(cfg.runtime)
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
    source_file = sac_path('src/vacusr.t.{0}'.format(cfg.usr_script))
    f = open(source_file, 'r')
    f_lines = f.readlines()

    for i,line in enumerate(f_lines):
        if line.strip() == "!### DELTA_X ###":
            f_lines[i+1] = ' delta_x = {0}d6\n'.format(cfg.delta_x)
        if line.strip() == "!### DELTA_Y ###":
            f_lines[i+1] = ' delta_y = {0}d6\n'.format(cfg.delta_y)
        if line.strip() == "!### DELTA_Z ###":
            f_lines[i+1] = ' delta_z = {0}d6\n'.format(cfg.delta_z)
        if line.strip() == "!### AMPLITUDE ###":
            f_lines[i+1] = '  AA = %s\n'%cfg.fort_amp
        if line.strip() == "!### PERIOD ###":
            f_lines[i+1] = '  s_period = %s\n'%cfg.period
        if line.strip() == "!### EXP_FAC ###":
            f_lines[i+1] = '  B = %s\n'%cfg.exp_fac
    f.close()

    #Truncate and overwrite
    f = open(source_file, 'w')
    f.writelines(f_lines)
    f.close()

    #==============================================================================
    # Compile Things
    #==============================================================================
    #Compile VAC
    os.chdir(sac_path("src"))
    os.system('./setvac -u=Slog -p=mhd -d=33 -g={0} -on=mpi'.format(cfg.grid_size))
    os.system('./setvac -s')
    if arguments['--clean']:
        os.system('./sac_fabricate.py --clean')
    else:
        os.system('./sac_fabricate.py')
    os.chdir(os.path.realpath('.'))
