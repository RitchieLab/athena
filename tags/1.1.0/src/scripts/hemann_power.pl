#!/usr/bin/perl

# updated GENN power script that looks for submodels of specified sizes
# within overall models

use strict;
$|=1;

my $PRGM_NAME = 'gennpower.pl';
my $VERSION_DATE = '1/21/09';

sub usage{
   print "\n\t$PRGM_NAME:  $VERSION_DATE\n\n";
   print "\tusage:\t\t$PRGM_NAME <loci> <max size> <outfile> <files>\n\n";
   print "\texample:\t$PRGM_NAME 5,10 3 genn.best.txt *.sum \n\n";
   print "\t<loci> are the loci of interest\n\t<max size> is the largest model size to consider\n\t<outfile> is the name of the file to store output (tab-delimited)\n\t<files> the hemann sum files to examine\n\n";
}

sub check_args{
  if(@ARGV < 1){
    usage();
    exit;
  }

  return @ARGV;
}

my @models;
my %sortmodels;
my $sortsize;

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


# returns array of models with each being the combo, total cv and
# average prediction error for the best in the indicated sizes
sub get_best{
  my ($file, $min, $max) = @_;
  
  @models= ();  
  
  open(IN, $file) or die "$PRGM_NAME: $file:  $!\n\n";
  while(<IN>){
    chomp;
    my @info = split" ";
    my $testval = $info[$#info];
    my @locs;
    for(my $i=1; $i < $#info - 1; $i++){
      # remove G from loc
      $info[$i] =~ s/G//;
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
    
    # generate combos for the sizes
    foreach my $size($min..$max){
      my @combos = combinations($size, scalar(@locs), 1);
      
      foreach my $combo(@combos){
        my @currcombo;
        # convert back to original loci values
        foreach my $val(@$combo){
          push(@currcombo, $locs[$val-1]);
        }
      
        my $key = join".",@currcombo;
        $models[$size]->{$key}->[0]++;
        $models[$size]->{$key}->[1] += $testval;
      }

    }    
  }
 
  close(IN);
  
  my @bestmodels;
  
  # select best models for each size
  foreach $sortsize($min..$max){
  
  # set error to be average error for all models
  foreach my $key(keys %{$models[$sortsize]}){
 #    print "$key cvc=$models[$sortsize]->{$key}->[0] pred=$models[$sortsize]->{$key}->[1]\n";
    $models[$sortsize]->{$key}->[1] /= $models[$sortsize]->{$key}->[0];
  }
    
    %sortmodels = %{$models[$sortsize]};
    my @sorted = sort by_cvc_testval keys %sortmodels;
    $bestmodels[$sortsize]->{key} = $sorted[0];
    $bestmodels[$sortsize]->{cvc} = $models[$sortsize]->{$sorted[0]}->[0];
    $bestmodels[$sortsize]->{pred} = $models[$sortsize]->{$sorted[0]}->[1];
  }
  
  return @bestmodels;
}

# 
sub by_cvc_testval{

  $sortmodels{$b}->[0] <=> $sortmodels{$a}->[0]
    ||
  $sortmodels{$a}->[1] <=> $sortmodels{$b}->[1]
}


my ($loci, $maxsize, $outfile, @files) = check_args();
my @loci = split/,/,$loci;

my $perfect_match=0; # power for perfect match
my $included_match=0; # total for match where model is included

# open each file and determine best model as selected
# first by best cross-validation consistency and then
# by testing error if there is a tie
my $total_files = @files;

sub best_by_cvc_pred{
  $b->{cvc} <=> $a->{cvc}
    ||
  $a->{pred} <=> $b->{pred}
}

open(OUT, ">$outfile") or die "$PRGM_NAME:  $!\n\n";

foreach my $size(2..$maxsize){
  print OUT "Best $size locus model\tCVC\tAvg Pred Error\t";
}

print OUT "Overall best model\tCVC\tAvg Pred Error or Test R2\n";


# counts number of times 
my @loci_found;

foreach my $file(@files){
  my @best = get_best($file, 2, $maxsize);
  
  my $bestcvc;
  my @bestcvc_size;
  
  # output each best model
  foreach my $size(2..$maxsize){
    my @locs = split/\./,$best[$size]->{key};
    print OUT "@locs\t";
    print OUT "$best[$size]->{cvc}\t";
    print OUT "$best[$size]->{pred}\t";
    
  }
  
  # show overall best
  my $bestmodel = (sort best_by_cvc_pred @best)[0];
    
  # output overall best model
  my @bestloci = split/\./,$bestmodel->{key};
  print OUT "@bestloci\t";
  print OUT "$bestmodel->{cvc}\t";
  print OUT "$bestmodel->{pred}\n";
  
  my %bestloci;
  my $numloci;
  
  foreach my $l(@bestloci){
    $bestloci{$l} = 1;
  }
  my $found;
  foreach my $loc(@loci){
    if($bestloci{$loc}){
      $found++;
      $numloci++;
    } 
  }
  #print "@bestloci\n";
  #print "found=$found\n";
  if($found == @loci){
    $included_match++;
    if($found == @bestloci){
      $perfect_match++;
    }
  }
  
  # add number found to each position
  for(my $i=$numloci; $i > 0; $i--){
    $loci_found[$i]++;
  }
  
}

print "\nMatches for @loci\n";
print "Number of files       = $total_files\n";
print "Perfect match         = $perfect_match\n";
print "Models including loci = $included_match\n";

my $total_num_loci = @loci;

for(my $i=$total_num_loci-1; $i > 0; $i--){
  print "Models containing at least $i of $total_num_loci loci   = $loci_found[$i]\n";
}
print "\n";
close(OUT);


