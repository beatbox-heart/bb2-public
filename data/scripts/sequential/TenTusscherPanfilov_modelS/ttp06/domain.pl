#!/usr/bin/perl -w
# The simulation domain as defined by fig 1(a), at given resolution
(($N=(($hx)=@ARGV))==1) or die "one arg please not $N\n";

$Lx=20.0;
$Ly=7.0;
$Lz=3.0;
$nx=$Lx/$hx;
$ny=$Ly/$hx;
$nz=$Lz/$hx;

for ($x=1;$x<=$nx+1;$x++) {
  for ($y=1;$y<=$ny+1;$y++) {
    for ($z=1;$z<=$nz+1;$z++) {
      print "$x, $y, $z, 1, 1, 0, 0\n";
    }
  }
}
