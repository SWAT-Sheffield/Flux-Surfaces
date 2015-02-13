#!/usr/bin/env python2
"""
Run code, analyis or build documents related to this repository

Usage:
    run.py SAC [--mpi --np=NP]
    run.py gdf [--mpi --np=NP]
    run.py fluxcalc [--mpi --np=NP]
    run.py surfflux (--tube-r=R) [--mpi --np=NP]
    run.py analysis (--tube-r=R) [--mpi --np=NP]
    run.py notebook [--port=PORT]
    run.py download (ini)

Options:
    --mpi  Use MPI, i.e. call mpiexec.
    --np=NP  Number of processors to use for MPI
    --tube-r=R  Specify the Flux Surface radius to use.
    --port=PORT Notebook port. [default: 8899]
"""
import os
import sys
import bz2
import subprocess

try:
    import docopt
except ImportError:
    from scripts.extern import docopt

try:
    import progressbar
except ImportError:
    from scripts.extern import progressbar

try:
    import requests
except ImportError:
    from scripts.extern import requests

os.chdir(os.path.dirname(os.path.realpath(__file__)))

arguments = docopt.docopt(__doc__, version='1.0.0')

#Fix rerun bug
if arguments['--rerun'] is None:
    arguments['--rerun'] = 'modified'

#MPI Exec function
def mpi_exec(arguments, command):
    if not arguments['--mpi']:
        os.system(command)
    if arguments['--mpi'] and arguments['--np']:
        os.system('mpirun -np %s %s'%(arguments['--np'], command))
    if arguments['--mpi'] and arguments['--np'] is None:
        os.system('mpirun %s'%command)

#Run SAC
if arguments['SAC'] == 'SAC':
    os.chdir("sac/sac")
    mpi_exec(arguments, './vac < vac.par')
    os.chdir("../")

#Run gdf translator
if arguments['gdf'] or arguments['SAC'] == 'gdf':
    os.chdir("sac/sac")
    mpi_exec(arguments, './out2gdf_pure.py')
    os.chdir("../")

#Run Analysis
if arguments['analysis']:
    os.chdir("analysis")
    mpi_exec(arguments, './surface_analysis_mpi.py --tube-r=%s'%arguments['--tube-r'])
    os.chdir("../")

#Run Flux Calcs
if arguments['surfflux']:
    os.chdir("analysis")
    mpi_exec(arguments, './surface_flux.py --tube-r=%s'%arguments['--tube-r'])
    os.chdir("../")

#Run Flux Calcs
if arguments['fluxcalc'] or arguments['SAC'] == 'fluxcalc':
    os.chdir("analysis")
    mpi_exec(arguments, './wave_flux.py')
    os.chdir("../")

#Start notebook sever
if arguments['notebook'] or arguments['SAC'] == 'notebook':
    if arguments['--port'] is None:
        arguments['--port'] = 8899
    import IPython
    IPython.start_ipython(['notebook', '--notebook-dir=analysis/notebooks',
    '--port={}'.format(arguments['--port'])])

#Download data
if arguments['download'] or arguments['SAC'] == 'download':
    ini_data_url = 'http://work.cadair.com/3D_tube_128_128_128.ini.bz'

    if arguments['ini']:
        data_url = ini_data_url
        filename = 'data/3D_tube_128_128_128.ini'
    else:
        sys.exit()

    CHUNK_SIZE = 1024 * 1024 # 1MB

    #Open file and get total size
    r = requests.get(data_url)
    total_size = int(r.headers['content-length'])
    #Create a progressbar
    print "Downloading and decompressing data:"
    pbar = progressbar.ProgressBar(widgets=[progressbar.Percentage(), progressbar.Bar()], maxval=total_size).start()
    #Create a stream bz2 decompressor
    bz = bz2.BZ2Decompressor()
    #Initial postition of file is 0
    pos = 0
    #Open the file and read and decompress chunk by chunk:
    with open(filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
            f.write(bz.decompress(chunk))
            #position is len of downloaded (compressed) file
            pos += len(chunk)
            pbar.update(pos)
    pbar.finish()
