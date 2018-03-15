#!/usr/bin/perl

use strict;

my $origfile = $ARGV[0] or die "$0 <datafile>";

my $continfile = $origfile . ".contin";

open(IN, $origfile) or die "$origfile:  $!\n\n";
open(OUT, ">$continfile") or die "$continfile:  $!\n\n";

srand(7);

while(<IN>){
  next unless /\d/;
  chomp;
  my @info = split" ";
  $info[0] = int(($info[0] * 20 + 50) * rand());   
  print OUT "@info\n";
}

close(OUT);
close(IN);

