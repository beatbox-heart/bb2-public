#!/bin/bash
# SK. August 2012.
# A bash shell script for preliminary sequential
# test runs of Beatbox in sequential mode (Beatbox_SEQ)
# using bbs scripts.
# 
rm -r out*
mkdir -p out
time $HOME/bin/Beatbox_SEQ LRD_3d_box.bbs
mv ppm ppm_LRD_3d_box
