#!/bin/bash
# SK. August 2012.
# A bash shell script for preliminary sequential
# test runs of Beatbox in sequential mode (bb2s)
# using bbs scripts.
# 
rm -r out*
mkdir -p out
time $HOME/bin/bb2s LRD_3d_box.bbs
mv ppm ppm_LRD_3d_box
