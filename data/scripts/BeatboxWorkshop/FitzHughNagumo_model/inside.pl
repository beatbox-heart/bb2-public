#!/usr/bin/perl -w
(($N=(($dir)=@ARGV))==1) or die "one argument please, not $N\n";

$task="$dir/task";
open TASK, ">$task" or die "cannot open $task for writing: $!";
print TASK "ppm
$dir/%04d.ppm
1,1,100                 mmin,mstep,mmax
123,113,135,3		nx,ny,nz (Number of gridpoints in each direction) 
1, 1 			normal resolution, rotation resolution

1                       show_sufrace
0, 0.35                 ulayer, uc: surface layer and constant
1, 0.6, 0.5             vlayer, vmin, vmax: red-colour layer and range
1, 0.6, 0.5             wlayer, wmin, wmax: blue-colour layer and range
3,0.1,1.0               color_mode, alphamin, alphamax

1                       show filament initially
0, 0.35                 filament layer 1 and const
1, 0.65                 filament layer 2 and const
1,1,0,1,2              	filament color, weight and balls flag

-75,70,0,5		theta, phi, psi, distance: initial view

0			back elimination initially
0			clipping plane initially

0                       write_images
0                       write_filaments
0,0,0,0                 write_history and the point

0,0,0                   background color
0.5,0.5,0.5,1           bounding box color and weight
1.0,1.0                 grid_h, dt: space and time quanta for dimensional output
3                       verbose: 0 for none, 1 for norm, 2 for detail, 3 for all

0			tracelen
0			sgzfile
0			tracez
0			traceint
0,0,0,0,0               color, alpha, weight most recent
0,0,0,0,0               color, alpha, weight least recent

0			markersize
0,0,0			marker coord
0,0,0,0,0		marker color alpha and weight

255			'taboo' byte value representing the void
";

system "ezview12 $task";
unlink $task;

