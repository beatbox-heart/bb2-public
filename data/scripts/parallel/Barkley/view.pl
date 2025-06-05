#!/usr/bin/perl -w
(($N=(($dir)=@ARGV))==1) or die "one argument please, not $N\n";

$task="$dir/task";
open TASK, ">$task" or die "cannot open $task for writing: $!";
print TASK "ppm
$dir/%06d.ppm
1,1,20             mmin,mstep,mmax
80,80,32,3		nx,ny,nz (Number of gridpoints in each direction) 
1, 1 			normal resolution, rotation resolution
0			field shown initially: 0 for u, 1 for v
0, 0.0, 1.0, 0.5     	ulayer,umin,umax,uc:     ! robust filament
1, 0.3, 0.8, 0.5        vlayer,vmin,vmax,vc:     ! detection 
1, 0.5, 0.8		color_mode, alphamin, alphamax
165, 60			theta, phi: initial view

0			back elimination initially
0			clipping plane initially
1			show filaments initially
0  			write filaments initially
0			save images initially
0  			write to history file initially
1,1,1                   (i,j,k) history point.   
1,1,1                    background colour
0.5,0.5,0.5             bounding box colour
0.2,0.8,0,1,1           filament colour, weight (and balls flag)
1.0,1.0			grid_h, dt: space and time quanta for dimensional output
0   			verbose: 0: none, 1: normal, 2: detailed, 3: all
50			markersize
100,100,34		marker coord
0,0,0,0,0		marker color alpha and weight
-1			'taboo' byte value representing the void
";

system "ezview11 $task";
unlink $task;

