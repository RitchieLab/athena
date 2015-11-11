#!/usr/bin/env perl

# Models will be selected as in MDR with only one model per CV.  Duplicate variables  will be 
# compressed into a single one.  So [G1 G1 G10 ]becomes [G1 G10]. 
# The best model from a set will be the one with highest CVC.  
# Balanced prediction accuracy will be the tie-breaker.
# Report the number of times the model of interest is the best one and the number of times
# it is either the best or included in the best model.  So if looking for [G1 G10] it will be
# included in [G1 G5 G10].
# For each smaller locus combinations that make up the model of interest, 
# report the number of times it is the best model and the number of times it appears
# in the best model.  So if [G1 G5 G10] is model of interest, report results on
# [G1 G5], [G1 G10], [G5, G10], [G1], [G5], [G10]


use strict;
use Getopt::Std;

$|=1;

my $PRGM_NAME = 'athena_power.pl';
my $VERSION_DATE = '12/29/11';

my %models;

sub usage{
  print "\n$PRGM_NAME:  $VERSION_DATE\n\n";
  print "usage:\t\t$PRGM_NAME [options] <variables> <files> \n\n";
  print "example:\t$PRGM_NAME G5,G10,C5 *.sum\n";
  print "\t\t$PRGM_NAME -l power.log G5,G10,C5 *.sum\n\n";
  print "options:\t-l\tLog file name listing all best models (default is no file)\n\n";
}


sub check_args{
  if(@ARGV < 2){
    usage();
    exit;
  }

  my %opts;
  getopts('l:', \%opts);

  my $logfile;
  if($opts{'l'}){
    $logfile = $opts{'l'};
  }

  return ($logfile, @ARGV);
}


sub by_cvc_testval{
  $models{$b}->{cvc} <=> $models{$a}->{cvc}
    ||
  $models{$b}->{pred} <=> $models{$a}->{pred}
}


# determines best model from the .sum file
sub get_best{
  my $file = shift;
  
  %models = (); 
  
  open(IN, $file) or die "$file:  $!\n\n";
  <IN>;
  while(<IN>){
    next unless /\d/;
    chomp;
    my @info = split" ";
    my $testval = $info[$#info];
    my @locs;
    for(my $i=1; $i < $#info - 1; $i++){
      # remove G from loc
#      $info[$i] =~ s/G//;
      push(@locs, $info[$i]);
    }
    
    my %unique_locs;
    # adapt locs back to original values
    foreach my $loc(@locs){
      $unique_locs{$loc}=1;
    }
    
    @locs = ();
    # remove any locs that are repeats
    foreach my $loc(sort keys %unique_locs){
      push(@locs, $loc);
    }
    my $key = join".",@locs;
    $models{$key}->{cvc}++;
    $models{$key}->{pred}+=$testval;
  }
  close(IN);
  
  #calculate prediction average and sort models
  foreach my $key(keys %models){
    $models{$key}->{pred} = $models{$key}->{pred}/$models{$key}->{cvc};
  }

  sub by_cvc_testval{
    $models{$b}->{cvc} <=> $models{$a}->{cvc}
      ||
    $models{$b}->{pred} <=> $models{$a}->{pred}
  }
 
  my @sorted = sort by_cvc_testval keys %models;

  #return best model
  return $sorted[0],$models{$sorted[0]}->{cvc},$models{$sorted[0]}->{pred};
}

# checks for perfect match and model within second model
sub check_match{
  my ($powermodel, $model) = @_;
  
  # contruct hash of model
  my %locs;
  foreach my $loc(split/\./, $model){
    $locs{$loc} = 1;
  }
  
  my @powerlocs = split/\./,$powermodel;
  
  my $matchcount=0;
  # count matches
  foreach my $powerloc(@powerlocs){
    if($locs{$powerloc}){
      $matchcount++;
    }
  }
  
  my $included=0;
  my $perfect_match = 0;
  
  # first check that all have been found
  if($matchcount==@powerlocs){
    $included = 1;
    if($matchcount==scalar(keys %locs)){
      $perfect_match = 1;
    }
  }
  return ($included, $perfect_match);
}

# recursive combination generator
# takes length of combination and number of elements
# returns the combinations
sub combinations{
  my $length = shift;
  my $n = shift;
  my $start = shift;
                    
  my @combos; 
  
  if($length == 0){
    return ();
  }
  if($length == 1){
    foreach my $i($start..$n){
      push(@combos, [$i]);
    }
    return @combos;
  }
  
  for my $i($start..$n){
    my @subcomb = combinations($length-1, $n, $i+1);
    foreach my $subcomb(@subcomb){
      unshift(@$subcomb, $i);
      push(@combos, $subcomb);
    }
  }
    
  return @combos;
}

my ($outfile, $loci, @files) = check_args();

my @loci = split/,/,$loci;

my $powermodel = join".", (sort @loci);

my $perfect_match=0; # power for perfect match
my $included_match=0; # total for match where model is included

my $total_files = @files;

if($outfile){
  open(OUT, ">$outfile") or die "$PRGM_NAME:  $!\n\n";
}
else{
  open(OUT, ">/dev/null")
}

print OUT "Best locus model\tCVC\tAvg Pred Score\n";

my @bestmodels;

my $bestmodelsize = 0;
my $nummods = 0;

foreach my $file(@files){
  my ($model, $modelcvc, $modelpred) = get_best($file);
  $nummods++;
  my @vars = split/\./, $model; 
  $bestmodelsize += scalar(@vars);

  push(@bestmodels, $model);
  my ($included, $perfect) = check_match($powermodel, $model);
  $included_match += $included;
  $perfect_match += $perfect;
    my @printlocs = split/\./,$model;
    print OUT "$file\t";
    print OUT "@printlocs\t";
    print OUT "$modelcvc\t";
    print OUT "$modelpred\n";
}
print OUT "\n\nAverage Best Model Size: ";
print OUT $bestmodelsize/$nummods;
print OUT "\n";
close(OUT);

print "\nMatches for @loci\n";
print "Number of files       = $total_files\n";
print "Perfect match         = $perfect_match\n";
print "Models including loci = $included_match\n\n";

my @pmodels;
# check each sub combination 
if(@loci > 1){
  foreach my $size(1..@loci-1){
    my @combos = combinations($size, scalar(@loci), 1);
    
    foreach my $combo(@combos){
      my @vars;
      foreach my $c(@$combo){
        push(@vars, $loci[$c-1]);
      }
      push(@pmodels, join '.', @vars);
    }
  }
}

# check each smaller model
for(my $i=$#pmodels; $i >=0; $i--){
  my $included_match=0;
  my $perfect_match=0;
  foreach my $bestmodel(@bestmodels){
    my ($included, $perfect) = check_match($pmodels[$i], $bestmodel);
    $included_match += $included;
    $perfect_match += $perfect;
  }
  
  my @printlocs = split/\./,$pmodels[$i];
  print "\nMatches for @printlocs\n";
  print "Number of files       = $total_files\n";
  print "Perfect match         = $perfect_match\n";
  print "Models including loci = $included_match\n\n";  
  
}


