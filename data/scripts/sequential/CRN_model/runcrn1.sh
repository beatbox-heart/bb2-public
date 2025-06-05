#!/bin/bash

# Sanjay Kharche.
# August 2012.
# Script for running the crn1.bbs script to get CV, CV restitution.

# the lowest and the highest premature pulse S2.
start=500
finish=600

for(( i = $start; i <= $finish; i=i+10 ))
do
S2=$i
fileName=crn1d$i.vtg
$HOME/bin/Beatbox_SEQ crn1.bbs $S2 $fileName 0.09 0.1235 
done
