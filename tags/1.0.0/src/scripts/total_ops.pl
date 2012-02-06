#!/usr/bin/perl

# totals up all appearances of the operators in 'best' files

use strict;
my $PRGM_NAME = 'total_ops.pl';


my %ops = ('PA'=>0, 'PD'=>0, 'PS'=>0, 'PM'=>0);

#subroutine: usage
# Displays current version date and usage
sub usage{
  print "\nusage:\t\t$PRGM_NAME <basefilename> \n\n";
  print "example:\t$PRGM_NAME results.200gen\n\n";
}

# subroutine: check_args
# Checks number of arguments passed in and returns name of pedigree file.
sub check_args{

  if( @ARGV < 1 )
  { usage();
    exit;}

  return @ARGV;
}

# parses best file from hemann and records operators used
sub parse_best{
  my $file = shift;
  
  open(IN, $file) or die "$file:  $!\n\n";
  while(<IN>){
    unless(/Grammar/){
        next;
    }
    $_ = <IN>;
    # next line is the network
    chomp;
    my @info = split" ";
    foreach my $terminal(@info){
      if(defined($ops{$terminal})){
        $ops{$terminal}++;
      }
    }
    last;
  }
  close(IN);
}


my ($basename) = check_args();

my @files = glob("$basename*.best");
foreach my $file(@files){
  parse_best($file);
}

print "Number of times operators appear in ";
print scalar @files;
print " files\n";
foreach my $op(keys %ops){
  print "$op $ops{$op}\n";
}

