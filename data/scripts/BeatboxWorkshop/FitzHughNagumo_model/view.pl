#!/usr/bin/perl -w
(($N=(($dir)=@ARGV))==1) or die "one argument please, not $N\n";

$task="$dir/task";
open TASK, ">$task" or die "cannot open $task for writing: $!";
print TASK "ppm
$dir/%06d.ppm
1,1,400             mmin,mstep,mmax

100,100,50,3		nx,ny,nz (Number of gridpoints in each direction) 
1, 1 			normal resolution, rotation resolution

1			show_surface
0,0.5			ulayer, uc
1,0.3,0.8		vlayer,vmin,vmax
1,0.3,0.8		wlayer,vmin,wmax
1, 0.5, 0.8		color_mode, alphamin, alphamax

1			show_filament
0,0.5			layer1,const1
1,0.5			layer2,const2
0.2,0.8,0,1,1           filament colour, weight (and balls flag)

165,60,0,5		theta,phi,psi,distance: initial view

0			back elimination initially
0			clipping plane initially

0			save images initially
0			save filament initially
0,1,1,1  		history flag and point

1,1,1                    background colour
0.5,0.5,0.5             bounding box colour
1.0,1.0			grid_h, dt: space and time quanta for dimensional output
0   			verbose: 0: none, 1: normal, 2: detailed, 3: all

0			tracelen
null			sgzfile
0                       tracez
0                       traceint
0,0,0,0,0,0,0,0		rgbaw,dx,dy,dz most recent
0,0,0,0,0,0,0,0		rgbaw,dx,dy,dz least recent

0			markersize
100,100,34		marker coord
0,0,0,0,0		marker color alpha and weight

-1			'taboo' byte value representing the void
";

system "ezview $task";
unlink $task;

