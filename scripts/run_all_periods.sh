#!/bin/bash

COMMAND='python ./scripts/swat_analysis.py'

./configure.py set driver --period=30.0 --amp=A20r2
$COMMAND
./configure.py set driver --period=60.0 --amp=A20
$COMMAND
./configure.py set driver --period=90.0 --amp=A20r2-3
$COMMAND
./configure.py set driver --period=120.0 --amp=A10R2
$COMMAND
./configure.py set driver --period=150.0 --amp=A4r10
$COMMAND
./configure.py set driver --period=180.0 --amp=A20-r3
$COMMAND
./configure.py set driver --period=210.0 --amp=A20r2-7
$COMMAND
./configure.py set driver --period=240.0 --amp=A10
$COMMAND
./configure.py set driver --period=270.0 --amp=A20-3r2
$COMMAND
./configure.py set driver --period=300.0 --amp=A4r5
$COMMAND
./configure.py set driver --period=330.0 --amp=A20r2-11
$COMMAND
