#!/usr/bin/env perl
# Make a .bbg geometry for the elliptic device tests.
#
# Usage: mkgeom.pl nx ny nz shape th0 th1 alpha [xcut] > file.bbg
#
#   nx,ny,nz  number of intervals in each direction (nz=0 for 2D, ny=nz=0 for 1D);
#             tissue nodes are x=1..nx+1 etc, so that x=0 and x=nx+2 are void.
#   shape     box   - the whole rectangular box is tissue
#             ball  - only nodes within the inscribed ellipse/ellipsoid are tissue
#   th0,th1   fibre angle in the xy plane, in degrees:
#             theta = th0 + th1*(x-1)/nx, i.e. th1 is the total rotation across the box
#   alpha     fibre elevation from the xy plane, degrees (constant)
#   xcut      optional: only the nodes with x<=xcut are tissue (the box or ball cut there)
#
# The fibre field must match the one assumed in ell_plane.bbs.
use strict;
use warnings;

die "Usage: $0 nx ny nz box|ball th0 th1 alpha [xcut] > file.bbg\n" if @ARGV < 7;
my ($nx, $ny, $nz, $shape, $th0, $th1, $alpha, $xcut) = @ARGV;
$xcut = $nx + 1 unless defined $xcut;
die "$0: shape must be box or ball\n" unless $shape eq 'box' || $shape eq 'ball';
my $deg = atan2(1, 1) / 45;

sub inside {
    my ($x, $y, $z) = @_;
    return 0 if $x > $xcut;
    return 1 if $shape eq 'box';
    my $r = 0.0;
    for my $cn ([$x, $nx], [$y, $ny], [$z, $nz]) {
        my ($c, $n) = @$cn;
        $r += (($c - 1 - $n / 2) / ($n / 2)) ** 2 if $n > 0;
    }
    return $r <= 1.0 + 1e-12;
}

my @ys = $ny > 0 ? (1 .. $ny + 1) : (0);
my @zs = $nz > 0 ? (1 .. $nz + 1) : (0);
for my $x (1 .. $nx + 1) {
    my $th = ($th0 + $th1 * ($x - 1) / $nx) * $deg;
    my $fx = cos($alpha * $deg) * cos($th);
    my $fy = cos($alpha * $deg) * sin($th);
    my $fz = sin($alpha * $deg);
    for my $y (@ys) {
        for my $z (@zs) {
            printf "%d, %d, %d, 1, %.17g, %.17g, %.17g\n", $x, $y, $z, $fx, $fy, $fz
                if inside($x, $y, $z);
        }
    }
}
