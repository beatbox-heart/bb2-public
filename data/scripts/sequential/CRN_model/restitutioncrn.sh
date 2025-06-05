#!/bin/bash
# Script for running the crn0.bbs script to get APD resitution of basal CRN human atrial cell model.
# This bash script loops over several values of the premature stimulus, after pulsing the cell model for 20 (Nstim)
# conditioning stimuli. Values of model parameters can also be passed to the script as shown by 
# gk1 and gcal. The output filenames are for time traces (fileName) and APD measurements (fileName2).
# APD is measured using the Poincare device as seen in the crn0_rest.bbs script.
# The output can be visualised using gnuplot (or any other plotting software). Contents
# of the .vtg files are visualised as (within gnuplot):
# plot 'crn_basal0.vtg' u 1:2 w l
# and those of the restitution as
# plot 'crn_basal.apd' u 1:2 w lp

start=100
finish=1000
incr=100

for(( i = $start; i <= $finish; i=i+$incr ))
do
Nstim=20
Period=$i
fileName=crn_basal$i.vtg
fileName2=crn_basal.apd
repolarisation=0.7
gk1=0.09
gcal=0.1235

$HOME/bin/Beatbox_SEQ crn0_rest.bbs $Nstim $Period $fileName $fileName2 $repolarisation $gk1 $gcal
done

