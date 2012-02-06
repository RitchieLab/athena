#!/usr/bin/perl

# add an id in ascending or descending order
# to the data files

use strict;

die "$0:  <data file> <asc or desc>" unless @ARGV == 2;

my $numlines = 0;
open(IN, $ARGV[0]) or die "$ARGV[0]:  $!\n\n";
while(<IN>){
  next unless /\d/;
  $numlines++;
}
close(IN);

open(OUT, ">$ARGV[0].id") or die "$!\n";

my $asc;

open(IN, $ARGV[0]) or die "$ARGV[0]:  $!\n\n";
my $id;
if($ARGV[1] =~ /asc/i){
  $asc = 1;
  $id = 1;
}
else{
  $asc = 0;
  $id = $numlines;
}

while(<IN>){
  next unless /\d/;
  chomp;
  my @info = split" ";
  unshift(@info, $id);
  if($asc){
    $id++;
  }
  else{
    $id--;
  }
  print OUT "@info";
  print OUT "\n";
}

close(OUT);

close(IN);