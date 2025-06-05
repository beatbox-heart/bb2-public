#!/usr/bin/perl -w

(($N=(($nx)=@ARGV))==1) or die "one arg please not $N\n";

$nx=~s/_.*$//;

$ny=$nx;

# xbl=1
# nx=xbr-xbl
# xbr=xbl+nx

for ($x=1;$x<=$nx+1;$x++) {
  for ($y=1;$y<=$ny+1;$y++) {
    print "$x, $y, 0, 1, 1, 0, 0\n";
  }
}

