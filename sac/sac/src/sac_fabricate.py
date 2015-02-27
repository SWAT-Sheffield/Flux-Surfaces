#! /usr/bin/env python
from fabricate import *
import sys
sys.path.append("../../../")
from scripts import sacconfig

cfg = sacconfig.SACConfig()

F_compiler = cfg.compiler
F_flags = cfg.compiler_flags.split(' ')
F_ext = '.f'

pre_processor = './vacpp.sh'
sources = ['vac', 'vacio', 'vacgrid', 'vacphys0', 'vacphys', 'vacusr']
includes = ['vacpar', 'vacusrpar', 'vacdef']
options = cfg.vac_modules

def build():
    pre_process()
    compile()
    link()

def pre_process():
    for source in sources + includes + options:
        run(pre_processor, source+'.t',  source+F_ext)

def compile():
    for source in sources + options:
        run(F_compiler, *F_flags + ['-c', source+F_ext])

def link():
    objects = [s+'.o' for s in sources + options]
    run(F_compiler, '-o', '../vac', objects)

def clean():
    autoclean()

if __name__ == '__main__':
    main()
