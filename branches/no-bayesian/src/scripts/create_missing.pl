#!/usr/bin/perl
use strict;

srand(7);

open(IN,$ARGV[0]) or die "$ARGV[0] $!\n\n";
my @lines = <IN>;
close(IN);

foreach my $line(@lines){
chomp($line);
my @info = split" ", $line;
for(my $i=1; $i<@info; $i++){
  if(rand() < .03){
    $info[$i] = -1;
  }
} 
my $newline = join" ",@info;
print "$newline\n";
}
