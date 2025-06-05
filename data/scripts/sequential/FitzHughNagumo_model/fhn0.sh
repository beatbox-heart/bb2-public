#!/bin/bash
# Script to simulate the FHN cell model AP.
cp $HOME/beatbox/data/parameters/fhn.par .
time $HOME/bin/Beatbox_SEQ fhn0.bbs
