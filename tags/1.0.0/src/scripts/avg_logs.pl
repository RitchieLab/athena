#!/usr/bin/perl

# compiles all the log files in a directory and
# outputs averages for each generation

use strict;

my $PRGM_NAME = 'avg_logs.pl';

my $GEN = 0;
my $VALIDNN = 1;
my $AVGSIZE = 2;
my $AVGFIT = 3;
my $MAXFIT = 4;
my $MINFIT = 5;
my $AVGGENO = 6;
my $AVGCOVAR = 7;
my $BESTMOD = 8;

# Subroutine: usage
# Displays current version date and usage
sub usage{
  print "\nusage:\t\t$PRGM_NAME <basefilename> \n\n";
  print "example:\t$PRGM_NAME results\n\n";
}

# subroutine: check_args
# Checks number of arguments passed in and returns name of pedigree file.
sub check_args{

  if( @ARGV < 1 )
  { usage();
    exit;}

  return @ARGV;
}


sub process_log{  
  my $name = shift;
  my $totals = shift;
  
  open(IN, $name) or die "$name:  $!\n\n";
  
  my $crossval;
  
  while(<IN>){
    if(/Cross\-validation\s#(\d+)/){
      $crossval = $1-1;
    }
    elsif(/^\d/){
      chomp;
      my @info = split" ";
      $totals->[$crossval]->[$info[$GEN]]->{'validnn'} += $info[$VALIDNN];
      $totals->[$crossval]->[$info[$GEN]]->{'avgsize'} += $info[$AVGSIZE];
      $totals->[$crossval]->[$info[$GEN]]->{'avgfit'} += $info[$AVGFIT];
      $totals->[$crossval]->[$info[$GEN]]->{'maxfit'} += $info[$MAXFIT];
      $totals->[$crossval]->[$info[$GEN]]->{'minfit'} += $info[$MINFIT];
      $totals->[$crossval]->[$info[$GEN]]->{'avggeno'} += $info[$AVGGENO];
      $totals->[$crossval]->[$info[$GEN]]->{'avgcovar'} += $info[$AVGCOVAR];
      
      my $gen = $info[$GEN];
      
      # remove everything after Best-mod
      @info = @info[$BESTMOD+1..$#info];
      
      # put together the two SNPs to return totals to original SNP
      for(my $i=0; $i<@info; $i++){
        $totals->[$crossval]->[$gen]->{'snptotals'}->[int($i/2)] += $info[$i];
      }
#  print "size of snptotals $gen = ";
#  print scalar @{$totals->[$crossval]->[$gen]->{'snptotals'}};
#  print "\n";
#  exit;
    }
  }
  close(IN);
#   last;
}


my ($basename) = check_args();
my $match = "$basename*.log";

my @files = glob($match);

my @totals;

foreach my $file(@files){
  process_log($file, \@totals);
}

my $totalfiles = scalar @files;

# output averages
for(my $cv=0; $cv < @totals; $cv++){
  my $outcv = $cv+1;
  print "Cross-validation #$outcv\n";
  print "Gen  Valid-NN   Avg-Size  Avg-Fit Max-Fit  Min-Fit  Avg-Geno Avg-Cov ";
  
  for(my $i=1; $i<= @{$totals[$cv][0]->{'snptotals'}}; $i++){
    print "$i ";
  }
  print "\n";
  
  
  for(my $gen=0; $gen < @{$totals[$cv]}; $gen++){
    printf "%-3d  %-7.1f    %-.2f     %-.3f   %-.3f    %-.3f    %-.2f     %-.2f ", $gen,
      $totals[$cv][$gen]->{'validnn'}/$totalfiles,
      $totals[$cv][$gen]->{'avgsize'}/$totalfiles, $totals[$cv][$gen]->{'avgfit'}/$totalfiles, 
      $totals[$cv][$gen]->{'maxfit'}/$totalfiles, $totals[$cv][$gen]->{'minfit'}/$totalfiles,
      $totals[$cv][$gen]->{'avggeno'}/$totalfiles, $totals[$cv][$gen]->{'avgcovar'}/$totalfiles;
      
      for(my $i=0; $i < @{$totals[$cv][$gen]->{'snptotals'}}; $i++){
        printf "%.2f ", $totals[$cv][$gen]->{'snptotals'}->[$i] / $totalfiles;
      }
    print "\n";
  }
  print "\n";
}




