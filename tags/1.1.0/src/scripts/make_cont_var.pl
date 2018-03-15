#!/usr/bin/perl

use strict;

# makes continuous variable values
# one is linked to the status 
# all others are random

my $TOTAL = 4;
my $LINKED = 2; # indexed of linked variable

if(@ARGV != 1){
  print "$0 <genotype file>\n";
  exit;
}

my @status;

open(IN, $ARGV[0]) or die "$ARGV[0]:  $!\n\n";
while(<IN>){
  next unless /\d/;
  my @info = split" ";
  push(@status, $info[0]); # store status in array
}
close(IN);

srand(1077);

open(OUT, ">$ARGV[0].contvar") or die "$!\n\n";

foreach my $st(@status){
  my @info;
  for(my $i=0; $i<$TOTAL; $i++){
    if($i!=$LINKED){
      $info[$i] = rand() * 100;
    }
    else{
      $info[$LINKED] = ($st * 30 + 40) * rand(); 
    }
  }
  
  print OUT "@info\n";
}

close(OUT);