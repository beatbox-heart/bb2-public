#!/usr/bin/perl

foreach $file (@ARGV) {
  $file=~/^(.*)_(.*)\.dat$/ or next;
  $hx=$1; $ht=$2;
  open FILE, $file;
  while (<FILE>) {
    @F=split;
    print "$hx , $ht & ".(join ' & ',@F)." \\\\\\hline\n";
  }
}
