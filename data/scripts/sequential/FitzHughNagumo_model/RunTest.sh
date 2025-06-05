#!/bin/bash
# A bash shell script for illustrative sequential test runs
# of Beatbox in sequential mode (Beatbox_SEQ)
# 
# 0D
time $HOME/bin/Beatbox_SEQ fhn0.bbs -profile
if which gnuplot>/dev/null; then
    $HOME/bin/Beatbox_SEQ fhn0_ap.bbs  -profile
    echo "plot 'fhn_ap.vtg' w l; pause 3 'you have 3 sec to view this ...'; quit; "| gnuplot
fi
# 1D
time $HOME/bin/Beatbox_SEQ fhn1.bbs -profile
time $HOME/bin/Beatbox_SEQ fhn1_load.bbs -profile
# 2D
time $HOME/bin/Beatbox_SEQ fhn2.bbs -profile
time $HOME/bin/Beatbox_SEQ singz.bbs -profile
time $HOME/bin/Beatbox_SEQ fhn_spiral_ffr_slice.bbs -profile
time $HOME/bin/Beatbox_SEQ fhn_spiral_ffr_slice_aniso.bbs -profile
time $HOME/bin/Beatbox_SEQ fhn_crossFieldStim_ffr_slice.bbs -profile
# 3D
time $HOME/bin/Beatbox_SEQ fhn3.bbs -profile
time $HOME/bin/Beatbox_SEQ fhn_ffr.bbs -profile


