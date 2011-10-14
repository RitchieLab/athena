#!/usr/bin/perl

use strict;

# reads through all log files and finds best model at each generation
# then produces output showing power for each generation of the run

my $PRGM_NAME = 'log_power.pl';

my $GEN = 0;
my $VALIDNN = 1;
my $AVGSIZE = 2;
my $AVGFIT = 3;
my $MAXFIT = 4;
my $MINFIT = 5;
my $AVGGENO = 6;
my $AVGCOVAR = 7;
my $BESTMOD = 8;

my %sortmodels;

# Subroutine: usage
# Displays current version date and usage
sub usage{
  print "\nusage:\t\t$PRGM_NAME <correct loci> <basefilename> \n\n";
  print "example:\t$PRGM_NAME 5,10 results\n\n";
}

# subroutine: check_args
# Checks number of arguments passed in and returns name of pedigree file.
sub check_args{

  if( @ARGV < 1 )
  { usage();
    exit;}

  return @ARGV;
}


# reads log file and returns best models for each generation
sub check_power{
  my $file = shift;
  my $powertotals = shift;
  my @loci = @{shift(@_)};
  
  my @models;
  my $cv;
  
  open(IN, $file) or die "$file:  $!\n\n";
  while(<IN>){
   if(/Cross\-validation\s#(\d+)/){
      $cv = $1-1;
    }
    elsif(/^\d/){
      chomp;
      my @info = split" ";
      $models[$info[$GEN]][$cv] = $info[$BESTMOD];
    }
  }
  close(IN);
  
  for(my $gen=0; $gen < @models; $gen++){
  # have all models now need to find best for each generation for each size
    my @best = select_best(2, 3, @{$models[$gen]});
    
    # show overall best
    my $bestmodel = (sort best_by_cvc_pred @best)[0];
    
    my @bestloci = split/\./,$bestmodel->{key};
# print "bestloci=> @bestloci loc => @loci\n";
    
    my %bestloci;
   foreach my $l(@bestloci){
     $bestloci{$l} = 1;
   }
   my $found;
   foreach my $loc(@loci){
     if($bestloci{$loc}){
      $found++;
     } 
   }
  #print "@bestloci\n";
  #print "found=$found\n";
    if($found == @loci){
      $powertotals->[$gen]->{'included_match'}++;
      if($found == @bestloci){
        $powertotals->[$gen]->{'perfect_match'}++;
      }
    }
# print "done with gen=$gen found= $found\n";
# exit;  
  }
  
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

sub select_best{
  my ($min, $max, @mods) = @_;
  
  my @models;
  
  foreach my $mod(@mods){
# print "model is $mod\n";
    my @locs = split/_/,$mod;
    
    foreach my $l(@locs){
      $l = int(($l+2)/2);  #convert back to original indexes
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

# print "combo is @locs\n";

    
    # generate combos for the sizes
    foreach my $size($min..$max){
      my @combos = combinations($size, scalar(@locs), 1);
      
      foreach my $combo(@combos){
        my @currcombo;
# print "combo is @$combo\n";
        # convert back to original loci values
        foreach my $val(@$combo){
          push(@currcombo, $locs[$val-1]);
        }
      
        my $key = join".",@currcombo;
# print "key=$key\n";
        $models[$size]->{$key}->[0]++;
        $models[$size]->{$key}->[1] += 0; # no prediction value on this log
      }

    }
  }
  
  my @bestmodels;

  
  # select best models for each size
  foreach my $sortsize($min..$max){
  
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


sub by_cvc_testval{

  $sortmodels{$b}->[0] <=> $sortmodels{$a}->[0]
    ||
  $sortmodels{$a}->[1] <=> $sortmodels{$b}->[1]
}


sub best_by_cvc_pred{
  $b->{cvc} <=> $a->{cvc}
    ||
  $a->{pred} <=> $b->{pred}
}


my ($loci, $basename) = check_args();


my @loci = split/,/,$loci;
my $match = "$basename*.log";

my @files = glob($match);
my $total_files = scalar @files;
my @powertotals;

foreach my $file(@files){
  check_power($file, \@powertotals, \@loci);
}

# print results for each generation
print "total runs=$total_files\n";
print "gen\tperfect\tincluded\n";
for(my $i=0; $i<@powertotals; $i++){
  print "$i\t$powertotals[$i]->{'perfect_match'}\t$powertotals[$i]->{'included_match'}\n";
}



