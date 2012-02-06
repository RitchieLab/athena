#!/usr/bin/perl

# compiles all the log files in a directory and
# outputs averages for each generation

use strict;

my $PRGM_NAME = 'all_cv_log.pl';

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
    elsif(/\d/){
      chomp;
      my @info = split" ";
      $totals->[$info[$GEN]]->{'validnn'} += $info[$VALIDNN];
      $totals->[$info[$GEN]]->{'avgsize'} += $info[$AVGSIZE];
      $totals->[$info[$GEN]]->{'avgfit'} += $info[$AVGFIT];
      $totals->[$info[$GEN]]->{'maxfit'} += $info[$MAXFIT];
      $totals->[$info[$GEN]]->{'minfit'} += $info[$MINFIT];
      $totals->[$info[$GEN]]->{'avggeno'} += $info[$AVGGENO];
      $totals->[$info[$GEN]]->{'avgcovar'} += $info[$AVGCOVAR];
      
      my $gen = $info[$GEN];
      
      # remove everything after Best-mod
      @info = @info[$BESTMOD+1..$#info];
      
      # put together the two SNPs to return totals to original SNP
      for(my $i=0; $i<@info; $i++){
        $totals->[$gen]->{'snptotals'}->[int($i/2)] += $info[$i];
      }
      
      
    }
  }
  close(IN);
  
}


my ($basename) = check_args();
my $match = "$basename*.log";

my @files = glob($match);

my @totals;

foreach my $file(@files){
  process_log($file, \@totals);
}

# assume 10 fold CV
my $totalfiles = scalar @files * 10;

# output averages
print "All Cross-validations combined\n";
print "Gen  Valid-NN   Avg-Size  Avg-Fit Max-Fit  Min-Fit  Avg-Geno Avg-Cov ";

for(my $i=1; $i<= @{$totals[0]->{'snptotals'}}; $i++){
    print "$i ";
 }
  print "\n";

for(my $gen=0; $gen < @totals; $gen++){
   printf "%-3d  %-7.1f    %-.2f     %-.3f   %-.3f    %-.3f    %-.2f     %-.2f ", $gen,
    $totals[$gen]->{'validnn'}/$totalfiles,
    $totals[$gen]->{'avgsize'}/$totalfiles, $totals[$gen]->{'avgfit'}/$totalfiles, 
    $totals[$gen]->{'maxfit'}/$totalfiles, $totals[$gen]->{'minfit'}/$totalfiles,
    $totals[$gen]->{'avggeno'}/$totalfiles, $totals[$gen]->{'avgcovar'}/$totalfiles;

    for(my $i=0; $i < @{$totals[$gen]->{'snptotals'}}; $i++){
      printf "%.2f ", $totals[$gen]->{'snptotals'}->[$i] / $totalfiles;
    }
    
   print "\n";  
    
}
  print "\n";





