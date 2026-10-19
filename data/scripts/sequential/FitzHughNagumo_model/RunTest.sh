#!/bin/bash
# A bash shell script for illustrative sequential test runs
# of Beatbox in sequential mode (bb2s)
# 
# 0D
time $HOME/bin/bb2s fhn0.bbs -profile
if which gnuplot>/dev/null; then
    $HOME/bin/bb2s fhn0_ap.bbs  -profile
    echo "plot 'fhn_ap.vtg' w l; pause 3 'you have 3 sec to view this ...'; quit; "| gnuplot
fi
# 1D
time $HOME/bin/bb2s fhn1.bbs -profile
time $HOME/bin/bb2s fhn1_load.bbs -profile
# 2D
time $HOME/bin/bb2s fhn2.bbs -profile
time $HOME/bin/bb2s singz.bbs -profile
time $HOME/bin/bb2s fhn_spiral_ffr_slice.bbs -profile
time $HOME/bin/bb2s fhn_spiral_ffr_slice_aniso.bbs -profile
time $HOME/bin/bb2s fhn_crossFieldStim_ffr_slice.bbs -profile
# 3D
time $HOME/bin/bb2s fhn3.bbs -profile
time $HOME/bin/bb2s fhn_ffr.bbs -profile


