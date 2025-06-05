#!/bin/bash
# SK. August 2012.
# A bash shell script for preliminary sequential
# test runs of Beatbox in sequential mode (Beatbox_SEQ)
# using bbs scripts.
# 
rm -r ppm*
mkdir -p ppm
cp $HOME/beatbox/data/parameters/fhn.par .
time $HOME/bin/Beatbox_SEQ fhn0.bbs
time $HOME/bin/Beatbox_SEQ fhn1.bbs
time $HOME/bin/Beatbox_SEQ fhn2.bbs
mv ppm ppm_fhn2
mkdir ppm
time $HOME/bin/Beatbox_SEQ fhn3.bbs
rm fhn.par
mv ppm ppm_fhn3
mkdir -p ppm
mkdir -p out
cp $HOME/beatbox/data/parameters/lrd.par .
time $HOME/bin/Beatbox_SEQ LRD_3d_box.bbs
mv ppm ppm_LRD_3d_box
